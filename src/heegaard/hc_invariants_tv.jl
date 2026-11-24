# cf: Regina 7.4.1, Triangulation3.turaevViro / Cyclotomic.field / __getitem__ — exact backend contract and rational coefficient transport.
import AkaghefAlgebra

struct ReginaTVError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO,e::ReginaTVError)=print(io,"ReginaTVError(",e.reason,"): ",e.detail)

"""Existing offline Regina 7.4.1 Python runner; no dependency installation.

On Windows this selects the already provisioned Ubuntu WSL consumer documented
in docs/ENVIRONMENT.md. Pass `runner::Cmd` to select another existing Python.
"""
function hc_regina_runner()
    if Sys.iswindows()
        wsl=Sys.which("wsl.exe")
        wsl===nothing && throw(ReginaTVError(:unavailable,"wsl.exe was not found"))
        return `$wsl -d Ubuntu --exec /home/akaghef/.local/bin/uv run --offline --no-project --with regina==7.4.1 python`
    end
    uv=Sys.which("uv")
    uv===nothing && throw(ReginaTVError(:unavailable,"uv was not found; supply an existing Python runner"))
    `$uv run --offline --no-project --with regina==7.4.1 python`
end

function _hctv_path(command::Cmd,path)
    if Sys.iswindows() && lowercase(basename(first(command.exec))) in ("wsl","wsl.exe")
        drive,tail=splitdrive(abspath(path))
        isempty(drive) && throw(ReginaTVError(:path,"WSL requires a drive-qualified script path"))
        return "/mnt/"*lowercase(drive[1:1])*replace(tail,'\\'=>'/')
    end
    path
end

function _hc_regina_tv_call(payload;runner=nothing,timeout::Real=120)
    isfinite(timeout) && timeout>0 || throw(ArgumentError("timeout must be finite and positive"))
    command=runner===nothing ? hc_regina_runner() : runner
    command isa Cmd || throw(ArgumentError("runner must be a Cmd ending in Python"))
    script=_hctv_path(command,joinpath(@__DIR__,"hc_invariants_tv.py"))
    mktempdir() do dir
        input,output,diagnostic=(joinpath(dir,n) for n in ("input.json","output.json","stderr.txt"))
        write(input,JSON3.write(payload))
        process=try
            backend_input=_hctv_path(command,input)
            run(pipeline(`$command $script --input $backend_input --timeout $timeout`,stdin=devnull,stdout=output,stderr=diagnostic);wait=false)
        catch e
            throw(ReginaTVError(:unavailable,sprint(showerror,e)))
        end
        if timedwait(()->process_exited(process),timeout+30;pollint=0.1)==:timed_out
            kill(process);wait(process)
            throw(ReginaTVError(:timeout,"runner exceeded $(timeout+30)s; native worker has its own $(timeout)s deadline; stderr: "*read(diagnostic,String)))
        end
        wait(process)
        result=try
            JSON3.read(read(output,String),Dict{String,Any})
        catch
            throw(ReginaTVError(:invalid_reply,read(diagnostic,String)))
        end
        success(process) || throw(ReginaTVError(:runner_failed,read(diagnostic,String)))
        if !get(result,"ok",false)
            e=get(result,"error",Dict{String,Any}())
            throw(ReginaTVError(Symbol(get(e,"type","backend_error")),string(get(e,"message",e))))
        end
        result
    end
end

"""Compute Regina's exact Turaev–Viro invariant from finite closed face pairings.

`r >= 3` (default 3) is Regina's parameter, corresponding to SU(2) level r-2.
For even r, `parity` is ignored and q₀ is the principal primitive 2r-th root.
For odd r, `parity=true` selects that same root; `false` selects the principal
primitive r-th root. Return an existing `AkaghefAlgebra.Cyclotomic` with its
explicit root specification; all coefficients are exact rationals.

The original TV normalization is retained: S³ = -(q₀-q₀⁻¹)²/(2r),
S²×S¹ = 1. No floating-point evaluation or further normalization is performed.
`simplify=true` uses Regina's topology-preserving triangulation simplifier.
`certificate=true` returns `(value, certificate)` including backend version,
field coefficients and tetrahedron counts. `timeout` bounds the native worker.
"""
function Z_TuraevViro(t::HCClosedTriangulation,r::Integer=3;
        parity::Bool=true,simplify::Bool=true,certificate::Bool=false,
        runner=nothing,timeout::Real=120)
    3<=r<=typemax(Int)÷2 || throw(ArgumentError("r must be an integer >= 3 with 2r fitting Int"))
    payload=hc_triangulation_data(t)
    merge!(payload,Dict("r"=>Int(r),"parity"=>parity,"simplify"=>simplify))
    answer=_hc_regina_tv_call(payload;runner=runner,timeout=timeout)
    value=_hctv_value(answer,r,parity)
    certificate ? (;value,certificate=answer) : value
end

function _hctv_value(answer,r,parity)
    c=Int(answer["field"])
    c==(isodd(r) && !parity ? r : 2r) || throw(ReginaTVError(:field,"unexpected conductor"))
    coefficients=map(answer["coefficients"]) do s
        parts=split(String(s),'/')
        length(parts) in (1,2) || throw(ReginaTVError(:coefficient,"invalid rational coefficient"))
        parse(BigInt,parts[1])//(length(parts)==1 ? big(1) : parse(BigInt,parts[2]))
    end
    raw=AkaghefAlgebra.TaftExact.Cyclotomic(c,coefficients)
    value=AkaghefAlgebra.Cyclotomic(raw,AkaghefAlgebra.RootSpec(c,1))
    value
end

"""Evaluate a nonempty vector of `(r, parity)` pairs in one native worker.

Return the corresponding vector of exact values, or `(values, certificate)`.
The worker deadline bounds the complete batch. This avoids repeated optional
backend startup when comparing several levels of the same triangulation.
"""
function Z_TuraevViro(t::HCClosedTriangulation,requests::AbstractVector;
        simplify::Bool=true,certificate::Bool=false,runner=nothing,timeout::Real=120)
    isempty(requests) && throw(ArgumentError("requests must be nonempty"))
    checked=Tuple{Int,Bool}[]
    for request in requests
        length(request)==2 || throw(ArgumentError("requests must contain (r, parity) pairs"))
        r,p=request
        r isa Integer && 3<=r<=typemax(Int)÷2 && p isa Bool || throw(ArgumentError("invalid (r, parity) pair"))
        push!(checked,(Int(r),p))
    end
    payload=hc_triangulation_data(t)
    merge!(payload,Dict("requests"=>checked,"simplify"=>simplify))
    answer=_hc_regina_tv_call(payload;runner=runner,timeout=timeout)
    results=answer["results"]
    length(results)==length(checked) || throw(ReginaTVError(:count,"incorrect batch result count"))
    values=[_hctv_value(a,r,p) for (a,(r,p)) in zip(results,checked)]
    certificate ? (;values,certificate=answer) : values
end

"""HC → existing checked forward triangulation → Regina exact Turaev–Viro.

Joining-disk selection and standard/general HC input follow `hc_triangulation`.
The selected disk is preserved in the certificate; conversion errors propagate.
"""
function Z_TuraevViro(hc::HeegaardCode,r::Integer=3;joining_disk=nothing,
        disk_name=nothing,disk_id=nothing,certificate::Bool=false,kwargs...)
    3<=r<=typemax(Int)÷2 || throw(ArgumentError("r must be an integer >= 3 with 2r fitting Int"))
    t=hc_triangulation(hc;joining_disk=joining_disk,disk_name=disk_name,disk_id=disk_id,signature=false)
    answer=Z_TuraevViro(t,r;certificate=certificate,kwargs...)
    if certificate
        answer.certificate["joining_disk"]=t.joining_disk
    end
    answer
end

function Z_TuraevViro(hc::HeegaardCode,requests::AbstractVector;joining_disk=nothing,
        disk_name=nothing,disk_id=nothing,certificate::Bool=false,kwargs...)
    t=hc_triangulation(hc;joining_disk=joining_disk,disk_name=disk_name,disk_id=disk_id,signature=false)
    answer=Z_TuraevViro(t,requests;certificate=certificate,kwargs...)
    certificate && (answer.certificate["joining_disk"]=t.joining_disk)
    answer
end

"""Evaluate several finite triangulations at the same vector of (r, parity) pairs.

One native worker handles the complete collection, under one `timeout`.
Return one vector of exact values per input, or `(values, certificate)`.
"""
function Z_TuraevViro(inputs::AbstractVector{<:HCClosedTriangulation},requests::AbstractVector;
        simplify::Bool=true,certificate::Bool=false,runner=nothing,timeout::Real=120)
    isempty(inputs) && throw(ArgumentError("inputs must be nonempty"))
    isempty(requests) && throw(ArgumentError("requests must be nonempty"))
    checked=Tuple{Int,Bool}[]
    for request in requests
        length(request)==2 || throw(ArgumentError("requests must contain (r, parity) pairs"))
        r,p=request
        r isa Integer && 3<=r<=typemax(Int)÷2 && p isa Bool || throw(ArgumentError("invalid (r, parity) pair"))
        push!(checked,(Int(r),p))
    end
    jobs=[merge(hc_triangulation_data(t),Dict("requests"=>checked,"simplify"=>simplify)) for t in inputs]
    answer=_hc_regina_tv_call(Dict("jobs"=>jobs);runner=runner,timeout=timeout)
    length(answer["jobs"])==length(inputs) || throw(ReginaTVError(:count,"incorrect input count"))
    values=map(answer["jobs"]) do job
        length(job["results"])==length(checked) || throw(ReginaTVError(:count,"incorrect request count"))
        [_hctv_value(a,r,p) for (a,(r,p)) in zip(job["results"],checked)]
    end
    certificate ? (;values,certificate=answer) : values
end

"""Bulk HC entry: preserve each input's own joining disk and exact forward path."""
function Z_TuraevViro(inputs::AbstractVector{<:HeegaardCode},requests::AbstractVector;
        certificate::Bool=false,kwargs...)
    tris=[hc_triangulation(hc;signature=false) for hc in inputs]
    answer=Z_TuraevViro(tris,requests;certificate=certificate,kwargs...)
    if certificate
        for (job,t) in zip(answer.certificate["jobs"],tris)
            job["joining_disk"]=t.joining_disk
        end
    end
    answer
end
