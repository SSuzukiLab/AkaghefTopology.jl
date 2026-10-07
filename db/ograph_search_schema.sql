PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS ograph_search_schema_meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT OR IGNORE INTO ograph_search_schema_meta(key, value)
VALUES ('schema_version', '1');
INSERT OR IGNORE INTO ograph_search_schema_meta(key, value)
VALUES ('state_contract', 'unweighted exact raw presentation');
INSERT OR IGNORE INTO ograph_search_schema_meta(key, value)
VALUES ('move_contract', 'MP, literal 02, and PS are distinct; H is replay-only');

CREATE TABLE IF NOT EXISTS ograph_search_states (
  state_id TEXT PRIMARY KEY REFERENCES states(state_id) ON DELETE CASCADE,
  state_version TEXT NOT NULL,
  raw_hash TEXT NOT NULL UNIQUE,
  geometry_hash TEXT NOT NULL,
  raw_gauss_json TEXT NOT NULL,
  orientation_json TEXT NOT NULL,
  real_gauss_json TEXT NOT NULL,
  real_orientation_json TEXT NOT NULL,
  real_vertex_count INTEGER NOT NULL CHECK (real_vertex_count >= 0)
);

CREATE TABLE IF NOT EXISTS ograph_search_runs (
  run_id TEXT PRIMARY KEY,
  script_path TEXT NOT NULL,
  script_sha256 TEXT NOT NULL,
  label TEXT NOT NULL,
  scope_json TEXT NOT NULL,
  status TEXT NOT NULL CHECK (
    status IN ('running', 'completed-found', 'completed-none', 'budget-exhausted', 'failed')
  ),
  counts_json TEXT NOT NULL,
  started_at_utc TEXT NOT NULL,
  finished_at_utc TEXT
);

CREATE TABLE IF NOT EXISTS ograph_search_expansions (
  expansion_id TEXT PRIMARY KEY,
  source_state_id TEXT NOT NULL REFERENCES ograph_search_states(state_id),
  move_kind TEXT NOT NULL CHECK (move_kind IN ('MP', '02', 'PS')),
  enumerator_version TEXT NOT NULL,
  selector_json TEXT NOT NULL,
  status TEXT NOT NULL CHECK (status IN ('building', 'completed', 'failed')),
  candidate_count INTEGER NOT NULL CHECK (candidate_count >= 0),
  success_count INTEGER NOT NULL CHECK (success_count >= 0),
  reject_count INTEGER NOT NULL CHECK (reject_count >= 0),
  error_count INTEGER NOT NULL CHECK (error_count >= 0),
  error_json TEXT,
  UNIQUE(source_state_id, move_kind, enumerator_version, selector_json)
);

CREATE TABLE IF NOT EXISTS ograph_search_run_expansions (
  run_id TEXT NOT NULL REFERENCES ograph_search_runs(run_id) ON DELETE CASCADE,
  expansion_id TEXT NOT NULL REFERENCES ograph_search_expansions(expansion_id),
  ordinal INTEGER NOT NULL CHECK (ordinal > 0),
  reused INTEGER NOT NULL CHECK (reused IN (0, 1)),
  PRIMARY KEY(run_id, expansion_id),
  UNIQUE(run_id, ordinal)
);

CREATE TABLE IF NOT EXISTS ograph_search_operations (
  operation_id TEXT PRIMARY KEY,
  source_state_id TEXT NOT NULL REFERENCES ograph_search_states(state_id),
  target_state_id TEXT REFERENCES ograph_search_states(state_id),
  move_kind TEXT NOT NULL CHECK (move_kind IN ('MP', '02', 'PS')),
  direction TEXT NOT NULL CHECK (direction IN ('forward', 'inverse')),
  engine_version TEXT NOT NULL,
  candidate_json TEXT NOT NULL,
  outcome TEXT NOT NULL CHECK (outcome IN ('success', 'rejected', 'error')),
  error_json TEXT,
  CHECK (
    (outcome = 'success' AND target_state_id IS NOT NULL AND error_json IS NULL)
    OR
    (outcome <> 'success' AND target_state_id IS NULL AND error_json IS NOT NULL)
  ),
  UNIQUE(source_state_id, engine_version, candidate_json)
);

CREATE TABLE IF NOT EXISTS ograph_search_expansion_operations (
  expansion_id TEXT NOT NULL REFERENCES ograph_search_expansions(expansion_id) ON DELETE CASCADE,
  ordinal INTEGER NOT NULL CHECK (ordinal > 0),
  operation_id TEXT NOT NULL REFERENCES ograph_search_operations(operation_id),
  PRIMARY KEY(expansion_id, ordinal),
  UNIQUE(expansion_id, operation_id)
);

-- Weight data is kept on replay trials, never on graph states or graph edges.
CREATE TABLE IF NOT EXISTS ograph_search_weighted_trials (
  trial_id TEXT PRIMARY KEY,
  run_id TEXT NOT NULL REFERENCES ograph_search_runs(run_id) ON DELETE CASCADE,
  source_state_id TEXT NOT NULL REFERENCES ograph_search_states(state_id),
  path_fingerprint TEXT NOT NULL,
  operation_ids_json TEXT NOT NULL,
  source_weights_json TEXT NOT NULL,
  h_moves_json TEXT NOT NULL,
  outcome TEXT NOT NULL CHECK (outcome IN ('admissible', 'rejected')),
  target_weights_json TEXT,
  failure_kind TEXT,
  failure_step INTEGER,
  message TEXT,
  observed_count INTEGER NOT NULL CHECK (observed_count > 0),
  CHECK (
    (outcome = 'admissible' AND target_weights_json IS NOT NULL
      AND failure_kind IS NULL AND failure_step IS NULL)
    OR
    (outcome = 'rejected' AND target_weights_json IS NULL
      AND failure_kind IS NOT NULL)
  ),
  UNIQUE(run_id, path_fingerprint, source_weights_json, h_moves_json)
);

CREATE INDEX IF NOT EXISTS ograph_search_states_geometry_idx
ON ograph_search_states(geometry_hash);

-- Exact feasibility / obstruction certificates after a geometric path is
-- selected. These do not introduce weighted states into the search graph.
CREATE TABLE IF NOT EXISTS ograph_search_constraint_trials (
  trial_id TEXT PRIMARY KEY,
  run_id TEXT NOT NULL REFERENCES ograph_search_runs(run_id),
  path_fingerprint TEXT NOT NULL,
  operation_ids_json TEXT NOT NULL,
  algorithm_version TEXT NOT NULL,
  result_json TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS ograph_search_constraint_trials_path_idx
ON ograph_search_constraint_trials(path_fingerprint,algorithm_version);

CREATE TABLE IF NOT EXISTS ograph_search_edge_weight_relations (
  operation_id TEXT NOT NULL REFERENCES ograph_search_operations(operation_id),
  algorithm_version TEXT NOT NULL,
  relation_json TEXT NOT NULL,
  PRIMARY KEY(operation_id,algorithm_version)
);

CREATE INDEX IF NOT EXISTS ograph_search_operations_source_idx
ON ograph_search_operations(source_state_id, move_kind, outcome);

CREATE INDEX IF NOT EXISTS ograph_search_operations_target_idx
ON ograph_search_operations(target_state_id);

CREATE INDEX IF NOT EXISTS ograph_search_weighted_trials_run_idx
ON ograph_search_weighted_trials(run_id, outcome);

CREATE VIEW IF NOT EXISTS successful_ograph_search_edges AS
SELECT operation_id, source_state_id, target_state_id, move_kind, direction,
       engine_version, candidate_json
FROM ograph_search_operations
WHERE outcome = 'success';
