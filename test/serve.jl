# Persistent Topology runner. Loads Revise and Topology once, then handles tests
# and script requests from test/client.sh. See test/TESTING.md for usage.
#
# Protocol (files in the queue directory owned by this server):
#   <id>.req   caller cwd followed by arguments, each terminated by NUL
#   <id>.log   the run's output
#   <id>.exit  the exit code: 0 pass, 1 fail, 3 Revise could not apply the edit
#   server.alive  touched every few seconds while the server waits
# When Revise cannot apply an edit (for example a changed struct), the server
# answers 3 before dispatching the request and exits.

const QUEUE = get(ENV, "TOPOLOGY_TEST_QUEUE", joinpath(@__DIR__, ".queue"))
mkpath(QUEUE)
const ALIVE = joinpath(QUEUE, "server.alive")
write(joinpath(QUEUE, "server.pid"), string(getpid()))

println("SERVE loading Revise and Topology"); flush(stdout)
using Revise
using Topology
println("SERVE preparing Revise baseline"); flush(stdout)
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

println("SERVE ready queue=$(QUEUE) pid=$(getpid())"); flush(stdout)
touch(joinpath(QUEUE, "server.ready"))

function run_script(cwd, path, args)
    script = normpath(joinpath(cwd, path))
    mod = Module(gensym(:TopologyExperiment))
    Core.eval(mod, :(const ARGS = $args))
    Core.eval(mod, :(const PROGRAM_FILE = $script))
    Core.eval(mod, :(include(p::AbstractString) = Base.include($mod, p)))
    Core.eval(mod, :(include(f::Function, p::AbstractString) = Base.include(f, $mod, p)))
    Core.eval(mod, :(eval(x) = Core.eval($mod, x)))
    cd(cwd) do
        Base.include(mod, script)
    end
    return 0
end

function reply(id, code)
    path = joinpath(QUEUE, id * ".exit")
    write(path * ".tmp", string(code))
    mv(path * ".tmp", path)
end

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
        payload = read(req, String)
        # Legacy, manually started servers accepted space-separated test arguments.
        cwd, args = if endswith(payload, '\0')
            fields = String.(split(payload, '\0'; keepempty=true))
            (fields[1], fields[2:end-1])
        else
            (pwd(), String.(split(strip(payload))))
        end
        rm(req; force=true)
        log = joinpath(QUEUE, id * ".log")
        code = 1
        @debug "Checking for source edits" request=id
        revised = try
            Revise.revise(throw=true); true
        catch err
            open(log, "w") do io
                println(io, "REVISE_FAILED: the server will exit before running the request.")
                showerror(io, err); println(io)
            end
            false
        end
        if !revised
            reply(id, 3)
            rm(ALIVE; force=true)
            exit(3)
        end
        @debug "Running request" request=id
        open(log, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    code = try
                        if !isempty(args) && first(args) == "--exec"
                            Base.invokelatest(run_script, cwd, args[2], args[3:end])
                        else
                            Base.invokelatest(main, args)
                        end
                    catch err
                        showerror(stdout, err, catch_backtrace()); println(); 1
                    end
                end
            end
        end
        reply(id, code)
    end
end

serve()
