#!/bin/sh
# Run test/run.jl or --exec <script.jl> [args...] on the persistent server.

T=$(CDPATH= cd "$(dirname "$0")" && pwd -P) || exit 1
P=$(dirname "$T")

native_pwd() {
  # Git Bash hands paths to native Windows Julia. Other shells use POSIX paths.
  pwd -W 2>/dev/null || pwd -P
}

if [ "$1" = "--serve-internal" ]; then
  lock=$2
  run=$3
  trap 'rm -f "$lock/queue"; rmdir "$lock"' 0
  TOPOLOGY_TEST_QUEUE=$(CDPATH= cd "$run" && native_pwd) || exit 1
  export TOPOLOGY_TEST_QUEUE
  cd "$P" || exit 1
  case $(uname -s) in
    MINGW*|MSYS*|CYGWIN*) set -- --pkgimages=no ;;
    *) set -- ;;
  esac
  "${TOPOLOGY_JULIA:-julia}" --startup-file=no "$@" --project=. test/serve.jl </dev/null > "$run/server.out" 2>&1
  code=$?
  printf '%s\n' "$code" > "$run/server.exit.tmp" && mv "$run/server.exit.tmp" "$run/server.exit"
  exit "$code"
fi

if [ "$1" = "--exec" ]; then
  if [ "$#" -lt 2 ] || [ ! -f "$2" ]; then
    echo 'Usage: client.sh --exec <existing script.jl> [args...]' >&2
    exit 2
  fi
  script_dir=$(CDPATH= cd "$(dirname -- "$2")" && native_pwd) || exit 1
  script=$script_dir/$(basename -- "$2")
  shift 2
  set -- --exec "$script" "$@"
fi
cwd=$(native_pwd) || exit 1
Q=${TOPOLOGY_TEST_QUEUE:-$T/.queue}
mkdir -p "$Q" || exit 1
Q=$(CDPATH= cd "$Q" && pwd -P) || exit 1
lock=$Q/server.lock
wait=${TOPOLOGY_TEST_WAIT:-3600}
case $wait in
  ''|*[!0-9]*) echo 'client: TOPOLOGY_TEST_WAIT must be a positive integer' >&2; exit 2 ;;
esac
[ "$wait" -gt 0 ] || exit 2
deadline=$(expr "$(date +%s)" + "$wait") || exit 2

wait_for_change() {
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "client: timed out after ${wait}s; $1" >&2
    exit 124
  fi
  sleep 1
}

request() {
  run=
  until [ -n "$run" ]; do
    if run=$(cat "$lock/queue" 2>/dev/null) && [ -n "$run" ]; then
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
      nohup sh "$T/client.sh" --serve-internal "$lock" "$run" </dev/null >/dev/null 2>&1 &
    else
      wait_for_change "waiting for server in $Q"
    fi
  done
  until [ -f "$run/server.ready" ] || [ -f "$run/server.exit" ]; do
    wait_for_change "waiting for server startup; see $run/server.out"
  done
  if [ -f "$run/server.exit" ]; then
    [ "$(cat "$run/server.exit")" = 3 ] && return 3
    cat "$run/server.out" >&2
    echo "client: Julia stopped before submission; see $run/server.out" >&2
    exit 1
  fi

  tmp=$(mktemp "$run/$(date +%Y%m%d%H%M%S)-request.XXXXXXXX") || exit 1
  printf '%s\000' "$cwd" "$@" > "$tmp" && mv "$tmp" "$tmp.req" || exit 1
  until [ -f "$tmp.exit" ]; do
    if [ -f "$run/server.exit" ]; then
      [ -f "$tmp.exit" ] && break
      if [ "$(cat "$run/server.exit")" = 3 ] && [ -f "$tmp.req" ]; then
        rm -f "$tmp.req"
        return 3
      fi
      [ ! -f "$tmp.log" ] || cat "$tmp.log"
      echo "client: Julia stopped before reporting a result; see $run/server.out" >&2
      exit 1
    fi
    wait_for_change "request may still run; see $tmp.log"
  done
  code=$(cat "$tmp.exit") || exit 1
  cat "$tmp.log" || exit 1
  rm -f "$tmp.log" "$tmp.exit"
  return "$code"
}

request "$@"
code=$?
if [ "$code" = 3 ]; then
  # Revise failed before dispatch, so this request can safely be submitted again.
  until [ "$(cat "$lock/queue" 2>/dev/null)" != "$run" ]; do
    wait_for_change "waiting for the old server to exit; see $run/server.out"
  done
  echo 'client: retrying on a fresh Julia server' >&2
  request "$@"
  code=$?
fi
exit "$code"
