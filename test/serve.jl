# Persistent test server. Loads Revise and Topology once, then runs test/run.jl
# requests from a queue directory, so a request pays neither the package load
# (about 180 s after a source edit) nor the first-call compilation of code it
# already ran. Clients use test/client.sh; policy in test/TESTING.md.
#
#   julia --project=. test/serve.jl        # from lib/topology.jl; client.sh starts it
#
# Protocol (files in the queue directory, default test/.queue):
#   <id>.req   one line of run.jl arguments, written by the client
#   <id>.log   the run's output
#   <id>.exit  the exit code: 0 pass, 1 fail, 3 Revise could not apply the edit
#   server.alive  touched every few seconds while the server waits
# When Revise cannot apply an edit (for example a changed struct), the server
# answers 3 and exits; the client then runs the request directly and the next
# client starts a fresh server.

using Revise
using Topology
includet(joinpath(@__DIR__, "run.jl"))

# Revise reads a package's source lazily from the precompile cache on the first edit.
# Lanes rebuild that cache whenever they load Topology after an edit, and then Revise has
# no baseline for the edited file and the server must exit. Parsing every file now keeps
# the baseline in memory, independent of the cache file.
let pkgdata = Revise.pkgdatas[Base.PkgId(Topology)], t = time()
    for file in Revise.srcfiles(pkgdata)
        Revise.maybe_parse_from_cache!(pkgdata, file)
    end
    println("SERVE parsed $(length(Revise.srcfiles(pkgdata))) Topology files in $(round(time() - t; digits=1)) s"); flush(stdout)
end

const QUEUE = get(ENV, "TOPOLOGY_TEST_QUEUE", joinpath(@__DIR__, ".queue"))
mkpath(QUEUE)
const ALIVE = joinpath(QUEUE, "server.alive")
write(joinpath(QUEUE, "server.pid"), string(getpid()))
println("SERVE ready queue=$(QUEUE) pid=$(getpid())"); flush(stdout)

function serve()
    last_beat = 0.0
    while true
        if time() - last_beat > 2
            touch(ALIVE); last_beat = time()
        end
        reqs = sort(filter(f -> endswith(f, ".req"), readdir(QUEUE)))
        if isempty(reqs)
            sleep(0.3); continue
        end
        req = joinpath(QUEUE, first(reqs))
        id = first(splitext(basename(req)))
        args = String.(split(strip(read(req, String))))
        rm(req; force=true)
        log = joinpath(QUEUE, id * ".log")
        code = 1
        revised = try
            Revise.revise(throw=true); true
        catch err
            open(log, "w") do io
                println(io, "REVISE_FAILED: the server cannot apply the current edit; run directly.")
                showerror(io, err)
            end
            false
        end
        if !revised
            write(joinpath(QUEUE, id * ".exit"), "3")
            rm(ALIVE; force=true)
            exit(3)
        end
        open(log, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    code = try
                        Base.invokelatest(main, args)
                    catch err
                        showerror(stdout, err, catch_backtrace()); 1
                    end
                end
            end
        end
        write(joinpath(QUEUE, id * ".exit"), string(code))
    end
end

serve()
