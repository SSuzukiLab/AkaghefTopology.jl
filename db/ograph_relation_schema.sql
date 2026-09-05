PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS ograph_schema_meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO ograph_schema_meta(key, value)
  VALUES ('schema_version', '2');

CREATE TABLE IF NOT EXISTS ograph_states (
  state_id TEXT PRIMARY KEY REFERENCES states(state_id) ON DELETE CASCADE,
  identity_domain TEXT NOT NULL,
  presentation_hash TEXT NOT NULL,
  convention_id TEXT NOT NULL,
  closedness_verifier_version TEXT NOT NULL,
  closedness_evidence TEXT NOT NULL,
  CHECK (length(trim(identity_domain)) > 0),
  CHECK (length(trim(closedness_evidence)) > 0)
);

CREATE TABLE IF NOT EXISTS ograph_sources (
  source_id TEXT PRIMARY KEY,
  source_kind TEXT NOT NULL,
  locator TEXT NOT NULL,
  snapshot_sha256 TEXT,
  metadata_json TEXT NOT NULL,
  UNIQUE(source_kind, locator, snapshot_sha256, metadata_json)
);

CREATE TABLE IF NOT EXISTS ograph_state_sources (
  state_id TEXT NOT NULL REFERENCES ograph_states(state_id) ON DELETE CASCADE,
  source_id TEXT NOT NULL REFERENCES ograph_sources(source_id) ON DELETE RESTRICT,
  PRIMARY KEY(state_id, source_id)
);

CREATE TABLE IF NOT EXISTS ograph_relation_profiles (
  profile_id TEXT PRIMARY KEY,
  profile_name TEXT NOT NULL CHECK (
    profile_name IN ('presentation', 'combed', 'manifold',
                     'integral-framed', 'framed-mod2', 'spin')
  ),
  profile_version INTEGER NOT NULL CHECK (profile_version > 0),
  parent_profile_id TEXT REFERENCES ograph_relation_profiles(profile_id),
  profile_fingerprint TEXT NOT NULL UNIQUE,
  moves_complete INTEGER NOT NULL CHECK (moves_complete IN (0, 1)),
  evidence_status TEXT NOT NULL CHECK (
    evidence_status IN ('identity-only', 'unimplemented', 'verified')
  ),
  conventions_json TEXT NOT NULL,
  created_at_utc TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS ograph_relation_profile_moves (
  profile_id TEXT NOT NULL
    REFERENCES ograph_relation_profiles(profile_id) ON DELETE CASCADE,
  move_type TEXT NOT NULL REFERENCES move_types(move_type),
  PRIMARY KEY(profile_id, move_type)
);

CREATE TABLE IF NOT EXISTS ograph_relation_campaigns (
  campaign_id TEXT PRIMARY KEY,
  profile_id TEXT NOT NULL
    REFERENCES ograph_relation_profiles(profile_id) ON DELETE RESTRICT,
  n_input INTEGER NOT NULL CHECK (n_input >= 0),
  b_path INTEGER NOT NULL CHECK (b_path >= n_input),
  input_count INTEGER NOT NULL CHECK (input_count >= 0),
  caps_json TEXT NOT NULL,
  input_fingerprint TEXT NOT NULL,
  campaign_fingerprint TEXT NOT NULL UNIQUE,
  engine_version TEXT NOT NULL,
  verifier_version TEXT NOT NULL,
  status TEXT NOT NULL CHECK (
    status IN ('planned', 'completed', 'budget-exhausted', 'failed')
  ),
  stop_reason TEXT NOT NULL,
  counts_json TEXT NOT NULL,
  created_at_utc TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS ograph_relation_campaign_inputs (
  campaign_id TEXT NOT NULL
    REFERENCES ograph_relation_campaigns(campaign_id) ON DELETE CASCADE,
  ordinal INTEGER NOT NULL CHECK (ordinal >= 0),
  state_id TEXT NOT NULL REFERENCES ograph_states(state_id) ON DELETE RESTRICT,
  PRIMARY KEY(campaign_id, ordinal),
  UNIQUE(campaign_id, state_id)
);

CREATE TABLE IF NOT EXISTS ograph_relation_certificates (
  certificate_id TEXT PRIMARY KEY,
  campaign_id TEXT NOT NULL
    REFERENCES ograph_relation_campaigns(campaign_id) ON DELETE RESTRICT,
  profile_id TEXT NOT NULL
    REFERENCES ograph_relation_profiles(profile_id) ON DELETE RESTRICT,
  source_state_id TEXT NOT NULL REFERENCES ograph_states(state_id),
  target_state_id TEXT NOT NULL REFERENCES ograph_states(state_id),
  certificate_fingerprint TEXT NOT NULL UNIQUE,
  status TEXT NOT NULL CHECK (status IN ('verified', 'rejected')),
  evidence_json TEXT NOT NULL,
  created_at_utc TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS ograph_relation_certificate_steps (
  certificate_id TEXT NOT NULL
    REFERENCES ograph_relation_certificates(certificate_id) ON DELETE CASCADE,
  ordinal INTEGER NOT NULL CHECK (ordinal >= 0),
  transition_id INTEGER NOT NULL REFERENCES transitions(id) ON DELETE RESTRICT,
  step_evidence_json TEXT NOT NULL,
  PRIMARY KEY(certificate_id, ordinal)
);

CREATE INDEX IF NOT EXISTS ograph_state_sources_source_idx
  ON ograph_state_sources(source_id);
CREATE INDEX IF NOT EXISTS ograph_campaigns_profile_idx
  ON ograph_relation_campaigns(profile_id);
CREATE INDEX IF NOT EXISTS ograph_certificates_profile_idx
  ON ograph_relation_certificates(profile_id);

CREATE TRIGGER IF NOT EXISTS ograph_states_domain_insert
BEFORE INSERT ON ograph_states
WHEN COALESCE((SELECT area FROM states WHERE state_id=NEW.state_id), -1.0) <> 3.0
  OR COALESCE((SELECT representation FROM states WHERE state_id=NEW.state_id), '') <> 'o-graph'
BEGIN SELECT RAISE(ABORT, 'o-graph states require Area 3 o-graph base states'); END;

CREATE TRIGGER IF NOT EXISTS ograph_states_domain_update
BEFORE UPDATE OF state_id ON ograph_states
WHEN COALESCE((SELECT area FROM states WHERE state_id=NEW.state_id), -1.0) <> 3.0
  OR COALESCE((SELECT representation FROM states WHERE state_id=NEW.state_id), '') <> 'o-graph'
BEGIN SELECT RAISE(ABORT, 'o-graph states require Area 3 o-graph base states'); END;

CREATE TRIGGER IF NOT EXISTS states_ograph_domain_update
BEFORE UPDATE OF area, representation ON states
WHEN EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.state_id)
 AND (NEW.area <> 3.0 OR NEW.representation <> 'o-graph')
BEGIN SELECT RAISE(ABORT, 'registered o-graph states must remain Area 3 o-graph states'); END;

CREATE TRIGGER IF NOT EXISTS ograph_profile_moves_area_insert
BEFORE INSERT ON ograph_relation_profile_moves
WHEN COALESCE((SELECT area FROM move_types WHERE move_type=NEW.move_type), -1.0) <> 3.0
BEGIN SELECT RAISE(ABORT, 'o-graph relation profiles require Area 3 moves'); END;

CREATE TRIGGER IF NOT EXISTS transitions_ograph_profile_insert
BEFORE INSERT ON transitions
WHEN (
   EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
   OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 )
 AND (
   NEW.relation_profile_id IS NULL
   OR NOT EXISTS (
     SELECT 1 FROM ograph_relation_profile_moves allowed
     WHERE allowed.profile_id=NEW.relation_profile_id
       AND allowed.move_type=NEW.move_type
   )
   OR NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
   OR NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 )
BEGIN SELECT RAISE(ABORT, 'Area 3 o-graph transitions require an allowed relation profile'); END;

CREATE TRIGGER IF NOT EXISTS transitions_ograph_profile_update
BEFORE UPDATE OF source_state_id, target_state_id, move_type, relation_profile_id ON transitions
WHEN (
   EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
   OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 )
 AND (
   NEW.relation_profile_id IS NULL
   OR NOT EXISTS (
     SELECT 1 FROM ograph_relation_profile_moves allowed
     WHERE allowed.profile_id=NEW.relation_profile_id
       AND allowed.move_type=NEW.move_type
   )
   OR NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
   OR NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 )
BEGIN SELECT RAISE(ABORT, 'Area 3 o-graph transitions require an allowed relation profile'); END;

CREATE TRIGGER IF NOT EXISTS transitions_nonograph_profile_insert
BEFORE INSERT ON transitions
WHEN NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
 AND NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 AND NEW.relation_profile_id IS NOT NULL
BEGIN SELECT RAISE(ABORT, 'non-o-graph transitions cannot carry an o-graph relation profile'); END;

CREATE TRIGGER IF NOT EXISTS transitions_nonograph_profile_update
BEFORE UPDATE OF source_state_id, target_state_id, relation_profile_id ON transitions
WHEN NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
 AND NOT EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 AND NEW.relation_profile_id IS NOT NULL
BEGIN SELECT RAISE(ABORT, 'non-o-graph transitions cannot carry an o-graph relation profile'); END;

-- Schema v2 has no trusted Area-3 transition writer.  Keep the generic table
-- physically closed until a later schema introduces a typed, replayable
-- MP/PS/CP verifier (including the nonlocal Method-B PS witness).
CREATE TRIGGER IF NOT EXISTS transitions_ograph_v2_disabled_insert
BEFORE INSERT ON transitions
WHEN NEW.relation_profile_id IS NOT NULL
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
BEGIN SELECT RAISE(ABORT, 'o-graph schema v2 does not accept relation transitions'); END;

CREATE TRIGGER IF NOT EXISTS transitions_ograph_v2_disabled_update
BEFORE UPDATE ON transitions
WHEN NEW.relation_profile_id IS NOT NULL
 OR OLD.relation_profile_id IS NOT NULL
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.target_state_id)
BEGIN SELECT RAISE(ABORT, 'o-graph schema v2 does not accept relation transitions'); END;

CREATE TRIGGER IF NOT EXISTS ograph_certificates_campaign_profile_insert
BEFORE INSERT ON ograph_relation_certificates
WHEN COALESCE((SELECT profile_id FROM ograph_relation_campaigns
               WHERE campaign_id=NEW.campaign_id), '') <> NEW.profile_id
BEGIN SELECT RAISE(ABORT, 'certificate and campaign relation profiles differ'); END;

CREATE TRIGGER IF NOT EXISTS ograph_campaign_inputs_bound_insert
BEFORE INSERT ON ograph_relation_campaign_inputs
WHEN COALESCE((SELECT complexity FROM states WHERE state_id=NEW.state_id), -1.0) >
     COALESCE((SELECT n_input FROM ograph_relation_campaigns
               WHERE campaign_id=NEW.campaign_id), -1)
BEGIN SELECT RAISE(ABORT, 'campaign input exceeds N_input'); END;

CREATE TRIGGER IF NOT EXISTS ograph_certificate_steps_profile_insert
BEFORE INSERT ON ograph_relation_certificate_steps
WHEN COALESCE((SELECT relation_profile_id FROM transitions
               WHERE id=NEW.transition_id), '') <>
     COALESCE((SELECT profile_id FROM ograph_relation_certificates
               WHERE certificate_id=NEW.certificate_id), '')
BEGIN SELECT RAISE(ABORT, 'certificate step relation profile differs'); END;

-- Area-3 evidence is append-only.  Corrections receive new content-derived
-- identifiers instead of rewriting or deleting the audit record in place.
CREATE TRIGGER IF NOT EXISTS ograph_states_immutable_update
BEFORE UPDATE ON ograph_states
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_states_immutable_delete
BEFORE DELETE ON ograph_states
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_sources_immutable_update
BEFORE UPDATE ON ograph_sources
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_sources_immutable_delete
BEFORE DELETE ON ograph_sources
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_state_sources_immutable_update
BEFORE UPDATE ON ograph_state_sources
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_state_sources_immutable_delete
BEFORE DELETE ON ograph_state_sources
BEGIN SELECT RAISE(ABORT, 'o-graph evidence rows are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_profiles_immutable_update
BEFORE UPDATE ON ograph_relation_profiles
BEGIN SELECT RAISE(ABORT, 'o-graph relation profiles are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_profiles_immutable_delete
BEFORE DELETE ON ograph_relation_profiles
BEGIN SELECT RAISE(ABORT, 'o-graph relation profiles are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_profile_moves_immutable_update
BEFORE UPDATE ON ograph_relation_profile_moves
BEGIN SELECT RAISE(ABORT, 'o-graph relation profiles are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_profile_moves_immutable_delete
BEFORE DELETE ON ograph_relation_profile_moves
BEGIN SELECT RAISE(ABORT, 'o-graph relation profiles are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_campaigns_immutable_update
BEFORE UPDATE ON ograph_relation_campaigns
BEGIN SELECT RAISE(ABORT, 'o-graph campaigns are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_campaigns_immutable_delete
BEFORE DELETE ON ograph_relation_campaigns
BEGIN SELECT RAISE(ABORT, 'o-graph campaigns are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_campaign_inputs_immutable_update
BEFORE UPDATE ON ograph_relation_campaign_inputs
BEGIN SELECT RAISE(ABORT, 'o-graph campaigns are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_campaign_inputs_immutable_delete
BEFORE DELETE ON ograph_relation_campaign_inputs
BEGIN SELECT RAISE(ABORT, 'o-graph campaigns are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_certificates_immutable_update
BEFORE UPDATE ON ograph_relation_certificates
BEGIN SELECT RAISE(ABORT, 'o-graph certificates are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_certificates_immutable_delete
BEFORE DELETE ON ograph_relation_certificates
BEGIN SELECT RAISE(ABORT, 'o-graph certificates are append-only'); END;

CREATE TRIGGER IF NOT EXISTS ograph_certificate_steps_immutable_update
BEFORE UPDATE ON ograph_relation_certificate_steps
BEGIN SELECT RAISE(ABORT, 'o-graph certificates are append-only'); END;
CREATE TRIGGER IF NOT EXISTS ograph_certificate_steps_immutable_delete
BEFORE DELETE ON ograph_relation_certificate_steps
BEGIN SELECT RAISE(ABORT, 'o-graph certificates are append-only'); END;

CREATE TRIGGER IF NOT EXISTS states_ograph_immutable_update
BEFORE UPDATE ON states
WHEN EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.state_id)
BEGIN SELECT RAISE(ABORT, 'registered o-graph base states are append-only'); END;
CREATE TRIGGER IF NOT EXISTS states_ograph_immutable_delete
BEFORE DELETE ON states
WHEN EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.state_id)
BEGIN SELECT RAISE(ABORT, 'registered o-graph base states are append-only'); END;

CREATE TRIGGER IF NOT EXISTS move_types_ograph_immutable_update
BEFORE UPDATE ON move_types
WHEN EXISTS (SELECT 1 FROM ograph_relation_profile_moves WHERE move_type=OLD.move_type)
BEGIN SELECT RAISE(ABORT, 'move types referenced by o-graph profiles are append-only'); END;
CREATE TRIGGER IF NOT EXISTS move_types_ograph_immutable_delete
BEFORE DELETE ON move_types
WHEN EXISTS (SELECT 1 FROM ograph_relation_profile_moves WHERE move_type=OLD.move_type)
BEGIN SELECT RAISE(ABORT, 'move types referenced by o-graph profiles are append-only'); END;

CREATE TRIGGER IF NOT EXISTS transitions_ograph_immutable_update
BEFORE UPDATE ON transitions
WHEN OLD.relation_profile_id IS NOT NULL
 OR NEW.relation_profile_id IS NOT NULL
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.target_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=NEW.target_state_id)
BEGIN SELECT RAISE(ABORT, 'o-graph transitions are append-only'); END;
CREATE TRIGGER IF NOT EXISTS transitions_ograph_immutable_delete
BEFORE DELETE ON transitions
WHEN OLD.relation_profile_id IS NOT NULL
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.source_state_id)
 OR EXISTS (SELECT 1 FROM ograph_states WHERE state_id=OLD.target_state_id)
BEGIN SELECT RAISE(ABORT, 'o-graph transitions are append-only'); END;
