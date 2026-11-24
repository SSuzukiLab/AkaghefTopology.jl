#!/bin/sh
# Run test/run.jl through the persistent test server (test/serve.jl).
#
#   sh lib/topology.jl/test/client.sh heegaard.koda
#   sh lib/topology.jl/test/client.sh --tier integration heegaard
#
# Arguments are those of test/run.jl. If no server is alive, this starts one in
# the background and runs the request directly this time. If the server cannot
# apply the current edit (exit 3), the request is also run directly.
# TOPOLOGY_TEST_WAIT caps the wait for an answer in seconds (default 3600).

T=$(cd "$(dirname "$0")" && pwd)          # lib/topology.jl/test
P=$(dirname "$T")                         # lib/topology.jl
Q=${TOPOLOGY_TEST_QUEUE:-$T/.queue}
mkdir -p "$Q"

alive() {  # server heartbeat newer than 10 s
  [ -f "$Q/server.alive" ] || return 1
  age=$(( $(date +%s) - $(date -r "$Q/server.alive" +%s) ))
  [ "$age" -lt 10 ]
}

direct() {
  msg=$1; shift
  echo "client: running directly ($msg)"
  cd "$P" && julia --project=. test/run.jl "$@"
}

if ! alive; then
  if [ ! -f "$Q/server.starting" ] || [ $(( $(date +%s) - $(date -r "$Q/server.starting" +%s) )) -gt 600 ]; then
    touch "$Q/server.starting"
    (cd "$P" && nohup julia --project=. test/serve.jl > "$Q/server.out" 2>&1 &)
  fi
  direct "no live server; one is starting for the next request" "$@"
  exit $?
fi
rm -f "$Q/server.starting"

id=$(date +%Y%m%d%H%M%S)-$$
printf '%s\n' "$*" > "$Q/$id.tmp" && mv "$Q/$id.tmp" "$Q/$id.req"
wait=${TOPOLOGY_TEST_WAIT:-3600}; t=0
while [ ! -f "$Q/$id.exit" ]; do
  sleep 1; t=$((t + 1))
  if [ "$t" -ge "$wait" ]; then echo "client: no answer within ${wait}s; request $id left in $Q"; exit 124; fi
  if [ $((t % 15)) -eq 0 ] && ! alive && [ ! -f "$Q/$id.exit" ]; then
    rm -f "$Q/$id.req"; direct "server stopped while waiting" "$@"; exit $?
  fi
done
code=$(cat "$Q/$id.exit")
cat "$Q/$id.log"
rm -f "$Q/$id.log" "$Q/$id.exit"
if [ "$code" = "3" ]; then direct "server could not apply the edit" "$@"; exit $?; fi
exit "$code"
