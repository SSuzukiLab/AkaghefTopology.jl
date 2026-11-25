#!/bin/sh
# Run test/run.jl or --exec <script.jl> [args...] on a persistent Julia server.
#
# Three tiers, chosen by the first argument (test/TESTING.md, "Tiers of the runner"):
#   --instant              near-instant requests; shared server; limit 15 s
#   (no tier argument)     work requests; shared server; limit 10 min
#   --long --limit <min>   a run of 10 minutes or more; its own Julia, which exits
#                          when the request ends. --batch is the old name for it
#                          (limit 60 min unless --limit is given).
# A request that passes the limit of its queue is stopped and its client exits
# 124. The shared server is then restarted for the requests queued behind it.

T=$(CDPATH= cd "$(dirname "$0")" && pwd -P) || exit 1
P=$(dirname "$T")

native_pwd() {
  # Git Bash hands paths to native Windows Julia. Other shells use POSIX paths.
  pwd -W 2>/dev/null || pwd -P
}

mtime() {
  # Prints nothing when the file is gone; callers must handle that.
  stat -c %Y "$1" 2>/dev/null || stat -f %m "$1" 2>/dev/null
}

is_windows() {
  case $(uname -s) in
    MINGW*|MSYS*|CYGWIN*) return 0 ;;
    *) return 1 ;;
  esac
}

stop_server_process() {
  # Stop exactly the Julia process that recorded its pid in the run directory.
  pid=$(tr -d '\r\n ' < "$1/server.pid" 2>/dev/null)
  case $pid in
    ''|*[!0-9]*) return 0 ;;
  esac
  if is_windows; then
    powershell -NoProfile -Command "\$p = Get-Process -Id $pid -ErrorAction SilentlyContinue; if (\$p -and \$p.ProcessName -eq 'julia') { Stop-Process -Id $pid -Force -Confirm:\$false }" >/dev/null 2>&1
  else
    kill "$pid" 2>/dev/null
  fi
}

if [ "$1" = "--serve-internal" ]; then
  # Supervisor of one server generation: starts Julia, enforces the limit of
  # the queue on every request, and releases the lock when Julia is gone.
  lock=$2
  run=$3
  cap=$4
  oneshot=$5
  hint=$6
  trap 'rm -f "$lock/queue"; rmdir "$lock" 2>/dev/null' 0
  touch "$run/server.super"
  printf '%s
' "$$" > "$run/supervisor.pid"
  TOPOLOGY_TEST_QUEUE=$(CDPATH= cd "$run" && native_pwd) || exit 1
  export TOPOLOGY_TEST_QUEUE
  if [ "$oneshot" = 1 ]; then
    TOPOLOGY_TEST_ONESHOT=1
    export TOPOLOGY_TEST_ONESHOT
  fi
  cd "$P" || exit 1
  if is_windows; then
    set -- --pkgimages=no
  else
    set --
  fi
  (
    "${TOPOLOGY_JULIA:-julia}" --startup-file=no "$@" --project=. test/serve.jl </dev/null > "$run/server.out" 2>&1
    printf '%s\n' "$?" > "$run/julia.code.tmp" && mv "$run/julia.code.tmp" "$run/julia.code"
  ) &
  until [ -f "$run/julia.code" ]; do
    touch "$run/server.super"
    for taken in "$run"/*.taken; do
      [ -f "$taken" ] || continue
      id=${taken%.taken}
      [ -f "$id.exit" ] && continue
      now=$(date +%s)
      message=
      # A marker can vanish between the test and the read when the request ends
      # at that moment; an unreadable time therefore means "not over the limit".
      if [ -f "$id.start" ]; then
        # The request is running: the limit of the queue applies.
        since=$(mtime "$id.start")
        [ -n "$since" ] || continue
        age=$((now - since))
        [ "$age" -gt "$cap" ] &&
          message="RUNNER: stopped after $age s; the limit of this queue is $cap s. $hint"
      else
        # The server is still applying source edits before the request.
        since=$(mtime "$taken")
        [ -n "$since" ] || continue
        age=$((now - since))
        [ "$age" -gt 300 ] &&
          message="RUNNER: stopped after $age s of applying source edits (limit 300 s), before the request ran."
      fi
      [ -f "$id.exit" ] && continue
      if [ -n "$message" ]; then
        printf '\n%s\n' "$message" >> "$id.log"
        # Stop Julia first: on Windows it holds the log open until it is gone.
        stop_server_process "$run"
        printf '124\n' > "$id.exit.tmp" && mv "$id.exit.tmp" "$id.exit"
        rm -f "$taken" "$id.start"
      fi
    done
    sleep 2
  done
  code=$(cat "$run/julia.code")
  printf '%s\n' "$code" > "$run/server.exit.tmp" && mv "$run/server.exit.tmp" "$run/server.exit"
  exit "$code"
fi

usage() {
  echo 'Usage: client.sh [--instant | --long --limit <minutes>] (--exec <script.jl> [args...] | <run.jl arguments>)' >&2
  exit 2
}

tier=work
limit=
legacy_batch=
while :; do
  case ${1:-} in
    --instant) tier=instant; shift ;;
    --long) tier=long; shift ;;
    --batch) tier=long; legacy_batch=1; shift ;;
    --limit)
      [ "$#" -ge 2 ] || usage
      limit=$2
      shift 2 ;;
    *) break ;;
  esac
done
case $limit in
  *[!0-9]*) echo 'client: --limit takes a whole number of minutes' >&2; exit 2 ;;
esac
case $tier in
  instant)
    cap=15
    default_wait=600
    hint='Send it without --instant.' ;;
  work)
    cap=600
    default_wait=3600
    hint='Shrink the case or split the question; start a run that must take longer with --long --limit <minutes>.' ;;
  long)
    if [ -z "$limit" ]; then
      [ -n "$legacy_batch" ] || { echo 'client: --long needs --limit <minutes>' >&2; exit 2; }
      limit=60
      echo 'client: --batch is now --long; using --limit 60' >&2
    fi
    [ "$limit" -gt 0 ] || { echo 'client: --limit must be positive' >&2; exit 2; }
    cap=$((limit * 60))
    default_wait=$((cap + 900))
    hint='Raise --limit only for an acceptance run.' ;;
esac
if [ "$tier" != long ] && [ -n "$limit" ]; then
  echo 'client: --limit belongs to --long' >&2
  exit 2
fi
if [ -n "${TOPOLOGY_TEST_CAP:-}" ]; then
  # Test hook: the limit, in seconds, of a server this call starts.
  case $TOPOLOGY_TEST_CAP in
    *[!0-9]*) echo 'client: TOPOLOGY_TEST_CAP must be a positive integer' >&2; exit 2 ;;
  esac
  cap=$TOPOLOGY_TEST_CAP
fi

if [ "$1" = "--exec" ]; then
  if [ "$#" -lt 2 ] || [ ! -f "$2" ]; then
    usage
  fi
  script_dir=$(CDPATH= cd "$(dirname -- "$2")" && native_pwd) || exit 1
  script=$script_dir/$(basename -- "$2")
  shift 2
  set -- --exec "$script" "$@"
fi
cwd=$(native_pwd) || exit 1
oneshot=0
if [ -n "${TOPOLOGY_TEST_QUEUE:-}" ]; then
  Q=$TOPOLOGY_TEST_QUEUE
  mkdir -p "$Q" || exit 1
elif [ "$tier" = long ]; then
  mkdir -p "$T/.queue/long" || exit 1
  Q=$(mktemp -d "$T/.queue/long/$(date +%Y%m%d%H%M%S).XXXXXX") || exit 1
else
  Q=$T/.queue/$tier
  mkdir -p "$Q" || exit 1
fi
[ "$tier" = long ] && oneshot=1
Q=$(CDPATH= cd "$Q" && pwd -P) || exit 1
lock=$Q/server.lock
wait=${TOPOLOGY_TEST_WAIT:-$default_wait}
case $wait in
  ''|*[!0-9]*) echo 'client: TOPOLOGY_TEST_WAIT must be a positive integer' >&2; exit 2 ;;
esac
[ "$wait" -gt 0 ] || exit 2
deadline=$(expr "$(date +%s)" + "$wait") || exit 2

wait_for_change() {
  now=$(date +%s)
  if [ "$now" -ge "$deadline" ]; then
    echo "client: timed out after ${wait}s; $1" >&2
    exit 124
  fi
  if [ "$1" != "${last_stage:-}" ]; then
    last_stage=$1
    stage_started=$now
    stage_reported=
  elif [ -z "$stage_reported" ] && [ "$((now - stage_started))" -ge 5 ]; then
    echo "client: $last_stage" >&2
    stage_reported=1
  fi
  sleep 1
}

supervisor_silent() {
  # The supervisor touches server.super every 2 s, also while Julia runs a request.
  [ -f "$1/server.exit" ] && return 1
  [ -d "$1" ] || return 0
  beat=$1/server.super
  [ -f "$beat" ] || beat=$1
  seen=$(mtime "$beat")
  [ -n "$seen" ] || return 1
  [ "$(( $(date +%s) - seen ))" -gt 20 ]
}

server_abandoned() {
  # True when the lock points at a server whose supervisor is gone. Checked
  # twice, so a machine waking from sleep is not mistaken for a dead server.
  supervisor_silent "$1" || return 1
  sleep 5
  supervisor_silent "$1"
}

discard_server() {
  echo "client: the server of $Q stopped without releasing its lock; discarding it" >&2
  stop_server_process "$1"
  [ -f "$1/server.exit" ] || printf '125\n' > "$1/server.exit"
  rm -f "$lock/queue"
  rmdir "$lock" 2>/dev/null
}

request() {
  run=
  started_here=
  until [ -n "$run" ]; do
    if run=$(cat "$lock/queue" 2>/dev/null) && [ -n "$run" ]; then
      if server_abandoned "$run"; then
        discard_server "$run"
        run=
        continue
      fi
      break
    fi
    if mkdir "$lock" 2>/dev/null; then
      if ! run=$(mktemp -d "$Q/run.XXXXXXXX") ||
         ! printf '%s\n' "$run" > "$lock/queue.tmp" ||
         ! mv "$lock/queue.tmp" "$lock/queue"; then
        rm -f "$lock/queue.tmp"
        rmdir "$lock"
        exit 1
      fi
      echo 'client: starting persistent Julia' >&2
      started_here=1
      nohup sh "$T/client.sh" --serve-internal "$lock" "$run" "$cap" "$oneshot" "$hint" </dev/null >/dev/null 2>&1 &
    else
      wait_for_change "waiting for server in $Q"
    fi
  done
  until [ -f "$run/server.ready" ] || [ -f "$run/server.exit" ]; do
    wait_for_change "waiting for server startup; see $run/server.out"
  done
  if [ -f "$run/server.exit" ]; then
    # A server another caller started has just stopped: submit to the next one.
    [ -z "$started_here" ] && return 3
    [ "$(cat "$run/server.exit")" = 3 ] && return 3
    cat "$run/server.out" >&2
    echo "client: Julia stopped before submission; see $run/server.out" >&2
    exit 1
  fi

  tmp=$(mktemp "$run/$(date +%Y%m%d%H%M%S)-request.XXXXXXXX") || exit 1
  printf '%s\000' "$cwd" "$@" > "$tmp" && mv "$tmp" "$tmp.req" || exit 1
  until [ -f "$tmp.exit" ]; do
    if [ ! -f "$run/server.exit" ] && supervisor_silent "$run" && server_abandoned "$run"; then
      discard_server "$run"
    fi
    if [ -f "$run/server.exit" ]; then
      [ -f "$tmp.exit" ] && break
      if [ -f "$tmp.req" ]; then
        # Never dispatched (Revise failure, or the request ahead was stopped):
        # it can safely be submitted again.
        rm -f "$tmp.req"
        return 3
      fi
      [ ! -f "$tmp.log" ] || cat "$tmp.log"
      echo "client: Julia stopped before reporting a result; see $run/server.out" >&2
      exit 1
    fi
    if [ -f "$tmp.req" ]; then
      wait_for_change "request queued in $Q; a submitted request may still run after the client times out"
    elif [ -f "$tmp.log" ]; then
      wait_for_change "running request (limit ${cap}s); see $tmp.log"
    else
      wait_for_change "request accepted; waiting for the source update check and execution setup"
    fi
  done
  code=$(cat "$tmp.exit") || exit 1
  cat "$tmp.log" || exit 1
  rm -f "$tmp.log" "$tmp.exit" 2>/dev/null
  return "$code"
}

attempt=1
request "$@"
code=$?
while [ "$code" = 3 ] && [ "$attempt" -lt 5 ]; do
  # Not dispatched, so the request can safely be submitted to the next server.
  until [ "$(cat "$lock/queue" 2>/dev/null)" != "$run" ]; do
    wait_for_change "waiting for the old server to exit; see $run/server.out"
  done
  echo 'client: retrying on a fresh Julia server' >&2
  attempt=$((attempt + 1))
  request "$@"
  code=$?
done
exit "$code"
