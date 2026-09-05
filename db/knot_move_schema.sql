PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS schema_meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO schema_meta(key, value) VALUES ('schema_version', '2');

CREATE TABLE IF NOT EXISTS states (
  id INTEGER PRIMARY KEY,
  state_id TEXT NOT NULL UNIQUE,
  state_hash TEXT NOT NULL,
  canonical_version INTEGER NOT NULL CHECK (canonical_version > 0),
  serialization TEXT NOT NULL,
  payload BLOB,
  area REAL NOT NULL CHECK (area IN (2.0, 2.5, 3.0)),
  representation TEXT NOT NULL,
  complexity REAL,
  UNIQUE(state_hash, canonical_version, serialization)
);

CREATE TABLE IF NOT EXISTS move_types (
  move_type TEXT PRIMARY KEY,
  area REAL NOT NULL CHECK (area IN (2.0, 2.5, 3.0)),
  inverse_type TEXT NOT NULL,
  inverse_direction TEXT NOT NULL CHECK (inverse_direction IN ('forward', 'inverse')),
  UNIQUE(move_type, area)
);

CREATE TABLE IF NOT EXISTS transitions (
  id INTEGER PRIMARY KEY,
  source_state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE,
  target_state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE,
  move_type TEXT NOT NULL REFERENCES move_types(move_type),
  direction TEXT NOT NULL CHECK (direction IN ('forward', 'inverse')),
  location TEXT,
  parameters TEXT,
  engine_version TEXT NOT NULL,
  verifier_version TEXT NOT NULL,
  verification_evidence TEXT,
  verified INTEGER NOT NULL CHECK (verified = 1),
  UNIQUE(source_state_id, target_state_id, move_type, direction, parameters)
);

CREATE TABLE IF NOT EXISTS certificate_paths (
  id INTEGER PRIMARY KEY,
  transition_id INTEGER NOT NULL REFERENCES transitions(id) ON DELETE CASCADE,
  path TEXT NOT NULL,
  UNIQUE(transition_id, path)
);

CREATE TABLE IF NOT EXISTS invariants (
  id INTEGER PRIMARY KEY,
  state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE,
  name TEXT NOT NULL,
  value TEXT NOT NULL,
  source TEXT,
  UNIQUE(state_id, name, value)
);

CREATE TABLE IF NOT EXISTS knotinfo_snapshots (
  id INTEGER PRIMARY KEY,
  snapshot_sha256 TEXT NOT NULL UNIQUE,
  source_url TEXT,
  schema_version INTEGER,
  raw_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS knotinfo_raw_records (
  snapshot_id INTEGER NOT NULL REFERENCES knotinfo_snapshots(id) ON DELETE CASCADE,
  source_row_index INTEGER NOT NULL,
  crossing_number INTEGER,
  name TEXT NOT NULL,
  headers_json TEXT NOT NULL,
  values_json TEXT NOT NULL,
  PRIMARY KEY(snapshot_id, source_row_index),
  UNIQUE(snapshot_id, name)
);

CREATE TABLE IF NOT EXISTS knot_types (
  name TEXT PRIMARY KEY,
  area REAL NOT NULL DEFAULT 2.0 CHECK (area = 2.0),
  state_id TEXT REFERENCES states(state_id),
  knotinfo_snapshot_id INTEGER,
  knotinfo_row_index INTEGER,
  invariants_json TEXT NOT NULL,
  FOREIGN KEY(knotinfo_snapshot_id, knotinfo_row_index)
    REFERENCES knotinfo_raw_records(snapshot_id, source_row_index)
);

CREATE INDEX IF NOT EXISTS transitions_source_idx ON transitions(source_state_id);
CREATE INDEX IF NOT EXISTS transitions_target_idx ON transitions(target_state_id);
CREATE INDEX IF NOT EXISTS invariants_state_idx ON invariants(state_id);

-- v2 classification metadata is dataset-level and intentionally one row per
-- method version.  State rows reference only that method version; no Carter
-- certificates, paths, or per-state polynomial values belong here.
CREATE TABLE IF NOT EXISTS classification_runs (
  method_version TEXT PRIMARY KEY,
  canonical_version INTEGER NOT NULL CHECK (canonical_version > 0),
  max_n INTEGER NOT NULL CHECK (max_n >= 0),
  barrier INTEGER NOT NULL CHECK (barrier >= max_n),
  census_fingerprint TEXT NOT NULL,
  classification_fingerprint TEXT NOT NULL UNIQUE,
  counts_json TEXT NOT NULL,
  conventions_json TEXT NOT NULL,
  created_at_utc TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS state_classifications (
  state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE,
  planar INTEGER NOT NULL CHECK (planar IN (0, 1)),
  knot_type TEXT,
  classification_status TEXT NOT NULL
    CHECK (classification_status IN ('classified', 'unresolved', 'nonclassical')),
  method_version TEXT NOT NULL
    REFERENCES classification_runs(method_version) ON DELETE CASCADE,
  PRIMARY KEY(method_version, state_id),
  CHECK (
    (planar = 0 AND knot_type IS NULL AND classification_status = 'nonclassical')
    OR
    (planar = 1 AND knot_type IS NULL AND classification_status = 'unresolved')
    OR
    (planar = 1 AND knot_type IS NOT NULL
      AND length(trim(knot_type)) > 0
      AND classification_status = 'classified')
  )
);

CREATE INDEX IF NOT EXISTS state_classifications_state_idx
  ON state_classifications(state_id);
CREATE INDEX IF NOT EXISTS state_classifications_method_idx
  ON state_classifications(method_version);

-- Classification is an Area-2 service.  Keep the Area-2/2.5/3 boundary
-- enforced in SQL as well as in the Julia API, including later state edits.
CREATE TRIGGER IF NOT EXISTS state_classifications_area_insert
BEFORE INSERT ON state_classifications
WHEN COALESCE((SELECT area FROM states WHERE state_id = NEW.state_id), -1.0) <> 2.0
BEGIN SELECT RAISE(ABORT, 'state classifications require Area 2 states'); END;

CREATE TRIGGER IF NOT EXISTS state_classifications_area_update
BEFORE UPDATE OF state_id ON state_classifications
WHEN COALESCE((SELECT area FROM states WHERE state_id = NEW.state_id), -1.0) <> 2.0
BEGIN SELECT RAISE(ABORT, 'state classifications require Area 2 states'); END;

CREATE TRIGGER IF NOT EXISTS states_classification_area_update
BEFORE UPDATE OF area ON states
WHEN NEW.area <> 2.0
 AND EXISTS (SELECT 1 FROM state_classifications WHERE state_id = OLD.state_id)
BEGIN SELECT RAISE(ABORT, 'classified states must remain in Area 2'); END;

CREATE TRIGGER IF NOT EXISTS transitions_verified_insert
BEFORE INSERT ON transitions WHEN NEW.verified <> 1
BEGIN SELECT RAISE(ABORT, 'only verified transitions may be persisted'); END;

CREATE TRIGGER IF NOT EXISTS transitions_evidence_insert
BEFORE INSERT ON transitions
WHEN NEW.verified = 1
 AND (NEW.verification_evidence IS NULL OR length(trim(NEW.verification_evidence)) = 0)
BEGIN SELECT RAISE(ABORT, 'verified transitions require verification evidence'); END;

CREATE TRIGGER IF NOT EXISTS transitions_evidence_update
BEFORE UPDATE OF verification_evidence, verified ON transitions
WHEN NEW.verified = 1
 AND (NEW.verification_evidence IS NULL OR length(trim(NEW.verification_evidence)) = 0)
BEGIN SELECT RAISE(ABORT, 'verified transitions require verification evidence'); END;

CREATE TRIGGER IF NOT EXISTS transitions_domain_insert
BEFORE INSERT ON transitions
WHEN (SELECT area FROM states WHERE state_id = NEW.source_state_id) <>
     (SELECT area FROM states WHERE state_id = NEW.target_state_id)
  OR (SELECT area FROM move_types WHERE move_type = NEW.move_type) <>
     (SELECT area FROM states WHERE state_id = NEW.source_state_id)
BEGIN SELECT RAISE(ABORT, 'transition source, target, and move areas must match'); END;

CREATE TRIGGER IF NOT EXISTS transitions_domain_update
BEFORE UPDATE OF source_state_id, target_state_id, move_type ON transitions
WHEN (SELECT area FROM states WHERE state_id = NEW.source_state_id) <>
     (SELECT area FROM states WHERE state_id = NEW.target_state_id)
  OR (SELECT area FROM move_types WHERE move_type = NEW.move_type) <>
     (SELECT area FROM states WHERE state_id = NEW.source_state_id)
BEGIN SELECT RAISE(ABORT, 'transition source, target, and move areas must match'); END;
