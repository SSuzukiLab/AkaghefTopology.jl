PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS hc_search_schema_meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO hc_search_schema_meta(key, value)
VALUES ('schema_version', '1');
INSERT OR IGNORE INTO hc_search_schema_meta(key, value)
VALUES ('state_contract', 'canonical HC form when the identity probe passes, otherwise the stored writing');
INSERT OR IGNORE INTO hc_search_schema_meta(key, value)
VALUES ('move_contract', 'witnesses from movable with omitted selections; twists use q=±1 on every diagram curve');

CREATE TABLE IF NOT EXISTS hc_search_states (
  state_id TEXT PRIMARY KEY REFERENCES states(state_id) ON DELETE CASCADE,
  state_version TEXT NOT NULL,
  identity_mode TEXT NOT NULL CHECK (identity_mode IN ('canonical', 'literal')),
  raw_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS hc_search_runs (
  run_id TEXT PRIMARY KEY,
  script_path TEXT NOT NULL,
  script_sha256 TEXT NOT NULL,
  label TEXT NOT NULL,
  scope_json TEXT NOT NULL,
  status TEXT NOT NULL CHECK (
    status IN ('running', 'completed', 'budget-exhausted', 'failed')
  ),
  counts_json TEXT NOT NULL,
  started_at_utc TEXT NOT NULL,
  finished_at_utc TEXT
);

CREATE TABLE IF NOT EXISTS hc_search_operations (
  operation_id TEXT PRIMARY KEY,
  run_id TEXT NOT NULL REFERENCES hc_search_runs(run_id) ON DELETE CASCADE,
  source_state_id TEXT NOT NULL REFERENCES hc_search_states(state_id),
  target_state_id TEXT REFERENCES hc_search_states(state_id),
  depth INTEGER NOT NULL CHECK (depth >= 0),
  row_id TEXT NOT NULL,
  engine_version TEXT NOT NULL,
  candidate_json TEXT NOT NULL,
  outcome TEXT NOT NULL CHECK (outcome IN ('success', 'rejected', 'error')),
  inverse_status TEXT NOT NULL CHECK (
    inverse_status IN ('ok', 'mismatch', 'error', 'unchecked')
  ),
  error_json TEXT,
  CHECK (
    (outcome = 'success' AND target_state_id IS NOT NULL)
    OR
    (outcome <> 'success' AND target_state_id IS NULL)
  ),
  UNIQUE(run_id, source_state_id, engine_version, candidate_json)
);

CREATE INDEX IF NOT EXISTS hc_search_operations_source_idx
ON hc_search_operations(source_state_id, row_id, outcome);

CREATE VIEW IF NOT EXISTS successful_hc_search_edges AS
SELECT operation_id, run_id, source_state_id, target_state_id, depth, row_id,
       engine_version, candidate_json, inverse_status
FROM hc_search_operations
WHERE outcome = 'success';
