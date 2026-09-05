PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS ograph_invariant_schema_meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO ograph_invariant_schema_meta(key, value)
  VALUES ('schema_version', '1');

CREATE TABLE IF NOT EXISTS ograph_pi1_presentations (
  presentation_id TEXT PRIMARY KEY,
  state_id TEXT NOT NULL REFERENCES ograph_states(state_id) ON DELETE RESTRICT,
  presentation_hash TEXT NOT NULL UNIQUE,
  algorithm_id TEXT NOT NULL,
  algorithm_version INTEGER NOT NULL CHECK (algorithm_version > 0),
  convention_id TEXT NOT NULL,
  generator_count INTEGER NOT NULL CHECK (generator_count >= 0),
  relator_count INTEGER NOT NULL CHECK (relator_count >= 0),
  relators_json TEXT NOT NULL CHECK (json_valid(relators_json)),
  presentation_serialization TEXT NOT NULL,
  evidence_json TEXT NOT NULL CHECK (json_valid(evidence_json)),
  created_at_utc TEXT NOT NULL,
  UNIQUE(state_id, algorithm_id, algorithm_version, convention_id),
  CHECK (length(presentation_hash) = 64),
  CHECK (length(trim(created_at_utc)) > 0)
);

CREATE TABLE IF NOT EXISTS ograph_h1_results (
  h1_result_id TEXT PRIMARY KEY,
  presentation_id TEXT NOT NULL
    REFERENCES ograph_pi1_presentations(presentation_id) ON DELETE RESTRICT,
  algorithm_id TEXT NOT NULL,
  algorithm_version INTEGER NOT NULL CHECK (algorithm_version > 0),
  backend TEXT NOT NULL,
  backend_version TEXT NOT NULL,
  free_rank INTEGER NOT NULL CHECK (free_rank >= 0),
  torsion_invariant_factors_json TEXT NOT NULL
    CHECK (json_valid(torsion_invariant_factors_json)),
  torsion_subgroup_order_decimal TEXT NOT NULL,
  h1_finite_order_decimal TEXT,
  h1_signature TEXT NOT NULL,
  result_fingerprint TEXT NOT NULL UNIQUE,
  evidence_json TEXT NOT NULL CHECK (json_valid(evidence_json)),
  created_at_utc TEXT NOT NULL,
  UNIQUE(presentation_id, algorithm_id, algorithm_version, backend, backend_version),
  CHECK (length(result_fingerprint) = 64),
  CHECK (length(torsion_subgroup_order_decimal) > 0
     AND torsion_subgroup_order_decimal NOT GLOB '*[^0-9]*'),
  CHECK (h1_finite_order_decimal IS NULL
     OR (length(h1_finite_order_decimal) > 0
     AND h1_finite_order_decimal NOT GLOB '*[^0-9]*')),
  CHECK ((free_rank = 0 AND h1_finite_order_decimal IS NOT NULL)
      OR (free_rank > 0 AND h1_finite_order_decimal IS NULL)),
  CHECK (length(trim(created_at_utc)) > 0)
);

-- Reserved optional enrichment.  Schema v1 has no public writer and the
-- replay verifier rejects nonempty rows rather than trusting backend labels.
CREATE TABLE IF NOT EXISTS ograph_group_backend_runs (
  backend_run_id TEXT PRIMARY KEY,
  presentation_id TEXT NOT NULL
    REFERENCES ograph_pi1_presentations(presentation_id) ON DELETE RESTRICT,
  backend TEXT NOT NULL,
  backend_version TEXT NOT NULL,
  operation TEXT NOT NULL,
  options_json TEXT NOT NULL CHECK (json_valid(options_json)),
  triviality_status TEXT NOT NULL CHECK (
    triviality_status IN ('PROVEN_TRIVIAL', 'PROVEN_NONTRIVIAL', 'UNKNOWN')
  ),
  finiteness_status TEXT NOT NULL CHECK (
    finiteness_status IN ('PROVEN_FINITE', 'PROVEN_INFINITE', 'UNKNOWN')
  ),
  pi1_order_decimal TEXT,
  external_group_id TEXT,
  result_json TEXT NOT NULL CHECK (json_valid(result_json)),
  evidence_json TEXT NOT NULL CHECK (json_valid(evidence_json)),
  transcript_sha256 TEXT,
  created_at_utc TEXT NOT NULL,
  CHECK (pi1_order_decimal IS NULL OR finiteness_status = 'PROVEN_FINITE'),
  CHECK (external_group_id IS NULL OR finiteness_status = 'PROVEN_FINITE'),
  CHECK (triviality_status <> 'PROVEN_TRIVIAL'
      OR (finiteness_status = 'PROVEN_FINITE' AND pi1_order_decimal = '1'))
);

CREATE TABLE IF NOT EXISTS ograph_finite_quotient_fingerprints (
  backend_run_id TEXT NOT NULL
    REFERENCES ograph_group_backend_runs(backend_run_id) ON DELETE RESTRICT,
  family TEXT NOT NULL,
  search_bound_json TEXT NOT NULL CHECK (json_valid(search_bound_json)),
  status TEXT NOT NULL CHECK (status IN ('FOUND', 'NOT_FOUND_WITHIN_BOUND', 'FAILED')),
  fingerprint_json TEXT NOT NULL CHECK (json_valid(fingerprint_json)),
  evidence_json TEXT NOT NULL CHECK (json_valid(evidence_json)),
  PRIMARY KEY(backend_run_id, family, search_bound_json)
);

CREATE TABLE IF NOT EXISTS ograph_manifold_identifications (
  identification_id TEXT PRIMARY KEY,
  state_id TEXT NOT NULL REFERENCES ograph_states(state_id) ON DELETE RESTRICT,
  backend_run_id TEXT NOT NULL
    REFERENCES ograph_group_backend_runs(backend_run_id) ON DELETE RESTRICT,
  manifold_identification TEXT NOT NULL CHECK (manifold_identification = 'S3'),
  identification_basis TEXT NOT NULL CHECK (identification_basis = 'trivial_pi1'),
  evidence_json TEXT NOT NULL CHECK (json_valid(evidence_json)),
  created_at_utc TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS ograph_pi1_state_idx
  ON ograph_pi1_presentations(state_id);
CREATE INDEX IF NOT EXISTS ograph_h1_bucket_idx
  ON ograph_h1_results(free_rank, h1_finite_order_decimal);
CREATE INDEX IF NOT EXISTS ograph_h1_signature_idx
  ON ograph_h1_results(h1_signature);

CREATE TRIGGER IF NOT EXISTS ograph_pi1_immutable_update
BEFORE UPDATE ON ograph_pi1_presentations
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_pi1_immutable_delete
BEFORE DELETE ON ograph_pi1_presentations
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_h1_immutable_update
BEFORE UPDATE ON ograph_h1_results
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_h1_immutable_delete
BEFORE DELETE ON ograph_h1_results
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_backend_runs_disabled_insert
BEFORE INSERT ON ograph_group_backend_runs
BEGIN SELECT RAISE(ABORT, 'external group backend replay is not implemented in schema v1'); END;
CREATE TRIGGER IF NOT EXISTS ograph_backend_runs_immutable_update
BEFORE UPDATE ON ograph_group_backend_runs
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_backend_runs_immutable_delete
BEFORE DELETE ON ograph_group_backend_runs
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_quotients_immutable_update
BEFORE UPDATE ON ograph_finite_quotient_fingerprints
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_quotients_immutable_delete
BEFORE DELETE ON ograph_finite_quotient_fingerprints
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_identifications_disabled_insert
BEFORE INSERT ON ograph_manifold_identifications
BEGIN SELECT RAISE(ABORT, 'S3 identification requires a future proof-replay backend'); END;
CREATE TRIGGER IF NOT EXISTS ograph_identifications_immutable_update
BEFORE UPDATE ON ograph_manifold_identifications
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_identifications_immutable_delete
BEFORE DELETE ON ograph_manifold_identifications
BEGIN SELECT RAISE(ABORT, 'o-graph derived invariants are append-only'); END;
