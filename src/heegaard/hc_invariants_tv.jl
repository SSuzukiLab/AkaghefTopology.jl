# cf: Regina 7.4.1, Triangulation3.turaevViro / Cyclotomic.field / __getitem__ — exact backend contract and rational coefficient transport.
import AkaghefAlgebra

"""
    ReginaTVError(reason, detail) -> ReginaTVError

Store a failure of the Regina Turaev–Viro adapter or its returned data.

# Arguments
`reason::Symbol` and `detail::String`, or convertible values. Conversion
errors propagate. No diagnostic is checked, and constructing the record
does not launch Regina or throw the exception.

# Returns
`ReginaTVError` with reason and explanatory text. `showerror` displays both.

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. The record is an adapter diagnostic, not a Turaev–Viro invariant value
   or a theorem distinguishing manifolds. [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
Construction does not verify that a backend failed. An error can concern
availability, execution or returned data; its existence does not imply that
the mathematical invariant is undefined or zero. The constructor installs
no dependency.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = ReginaTVError(:backend, "backend unavailable"); (e.reason, e.detail)
(:backend, "backend unavailable")
```

See also [`Z_TuraevViro`](@ref), [`hc_regina_runner`](@ref).

Role: construct
"""
struct ReginaTVError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO,e::ReginaTVError)=print(io,"ReginaTVError(",e.reason,"): ",e.detail)

"""
    hc_regina_runner() -> Cmd

Return the command selecting the repository's offline Regina Python environment.

# Arguments
The Regina environment documented in `docs/ENVIRONMENT.md` must
already be provisioned. Checked: `wsl.exe` on Windows or `uv` elsewhere is
on `PATH`; absence raises `ReginaTVError(:unavailable)` explaining the
missing launcher. Backend availability itself is not checked.

# Returns
A `Cmd` suitable for the Regina adapter's `runner` keyword.

# Guarantees
1. The command requests the pinned Regina 7.4.1 environment in offline mode.
   [proved/agent: by construction]
2. The function launches no process and modifies no input or environment.
   [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
Returning a command does not verify that its Python environment, WSL distribution
or Regina package exists. Later command execution can fail; this function
computes no topological invariant.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — external-backend boundary; machine provisioning is recorded in docs/ENVIRONMENT.md.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.hc_regina_runner() isa Cmd catch e; e isa H.ReginaTVError && e.reason == :unavailable end
true
```

Role: other
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

"""
    Z_TuraevViro(t::HCClosedTriangulation, r::Integer=3; parity=true, simplify=true, certificate=false, runner=nothing, timeout=120) -> Union{AkaghefAlgebra.Cyclotomic,NamedTuple}

Return Regina's exact Turaev–Viro value for a supplied finite closed triangulation.

# Arguments
`t` supplies the finite pairings of [the triangulation interface](../specs/heegaard-snappy.md).
Regina checks the backend's accepted triangulation domain. Checked locally:
`3 <= r <= typemax(Int)÷2`, finite positive `timeout`, and `runner` either
`nothing` or a command ending in an existing Python interpreter; violations
raise `ArgumentError`. `parity`, `simplify` and `certificate` are Boolean.
Unavailable backend, timeout, rejected input and malformed replies raise
`ReginaTVError` identifying the failure. No dependency is installed here.

# Returns
An exact cyclotomic value, or `(value, certificate)` when requested. The
certificate contains the backend reply, including version, field coefficients
and triangulation counts.

# Guarantees
1. The result preserves the exact rational coefficients and root field
   supplied by the backend for this request. [proved/agent: by construction]
2. For an accepted triangulation, the intended invariant is that of its closed
   3-manifold in Regina's requested parameter and parity convention.
   [unverified: general source alignment of the adapter and normalization]
3. The input triangulation record is not modified. [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
The returned certificate does not independently identify the triangulation
with an HC's manifold. The [normalization specification](../specs/heegaard-code.md#moo-tv-normalization)
states the root and parity correspondence. Runtime depends on the supplied
triangulation; the process watchdog includes backend-cleanup grace time.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> t = HCClosedTriangulation(0, NamedTuple[], Int[], nothing, nothing, Dict{String,Any}(), NamedTuple()); try Z_TuraevViro(t, 2) catch e; e isa ArgumentError end
true
```

See also [`Z_TuraevViro_skeleton`](@ref), [`hc_regina_runner`](@ref).

Role: invariant
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

"""
    Z_TuraevViro(t::HCClosedTriangulation, requests::AbstractVector; simplify=true, certificate=false, runner=nothing, timeout=120) -> Union{Vector,NamedTuple}

Return exact Turaev–Viro values for several parameters on one finite triangulation.

# Arguments
`t`, `simplify`, `runner` and `timeout` follow the scalar triangulation overload.
`requests` must be nonempty, with each entry an `(r, parity)` pair, an integer
`r` in the scalar range and a Boolean parity. Checked violations raise
`ArgumentError` naming an empty request list or malformed pair. Backend and
conversion errors propagate as in the scalar overload; a wrong reply count
raises `ReginaTVError(:count)`. `certificate` is Boolean.

# Returns
A vector in request order, or `(values, certificate)` containing that vector
and the complete backend reply.

# Guarantees
1. Each value retains the exact backend coefficients for the corresponding
   request on the supplied triangulation. [proved/agent: by construction]
2. The input and request vector are not modified. [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
The scalar overload's geometric-identification and normalization limits
apply. The timeout covers the complete batch, not each request separately.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> t = HCClosedTriangulation(0, NamedTuple[], Int[], nothing, nothing, Dict{String,Any}(), NamedTuple()); try Z_TuraevViro(t, Tuple{Int,Bool}[]) catch e; e isa ArgumentError end
true
```

Role: invariant
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

"""
    Z_TuraevViro(hc::HeegaardCode, r::Integer=3; joining_disk=nothing, disk_name=nothing, disk_id=nothing, certificate=false, kwargs...) -> Union{AkaghefAlgebra.Cyclotomic,NamedTuple}

Return the exact Turaev–Viro value of the forward triangulation of an HC.

# Arguments
`hc` and disk selectors must satisfy `hc_triangulation`; its conversion
failures propagate. The scalar range for `r` is checked before conversion,
with `ArgumentError` on failure. Remaining keywords follow the scalar
triangulation overload; unsupported keywords raise `MethodError`.

# Returns
An exact cyclotomic value, or `(value, certificate)`. The requested certificate
also records the forward triangulation's selected joining disk.

# Guarantees
1. The result evaluates the triangulation produced by the existing HC
   forward conversion in [the triangulation interface](../specs/heegaard-snappy.md).
   [proved/agent: by construction]
2. The HC is not modified. [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
Identification with the intended HC 3-manifold inherits the forward
conversion's evidence limits. The scalar triangulation overload's backend
and normalization limitations apply.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(0, Dict(), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc, Dict(), Dict(), Dict(), Dict(), nothing); try Z_TuraevViro(h, 2) catch e; e isa ArgumentError end
true
```

See also [`hc_triangulation`](@ref).

Role: invariant
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

"""
    Z_TuraevViro(hc::HeegaardCode, requests::AbstractVector; joining_disk=nothing, disk_name=nothing, disk_id=nothing, certificate=false, kwargs...) -> Union{Vector,NamedTuple}

Return several exact Turaev–Viro values for the forward triangulation of an HC.

# Arguments
`hc` and disk selectors follow `hc_triangulation`; its errors propagate.
`requests`, `certificate` and remaining keywords follow the single-triangulation
batch overload. That overload checks the request pairs and backend reply.
Conversion occurs before request validation, so a bad HC can fail first.

# Returns
A request-ordered vector, or `(values, certificate)` whose certificate also
records the selected joining disk.

# Guarantees
1. Every value evaluates the same forward triangulation at its corresponding
   request. [proved/agent: by construction]
2. The HC and request vector are not modified. [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
The forward conversion's manifold-identification limits and the triangulation
overload's normalization and backend limits apply to this wrapper.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(-1, Dict(), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc, Dict(), Dict(), Dict(), Dict(), nothing); try Z_TuraevViro(h, [(3,true)]) catch e; e isa KodaInputError && e.stage == :validate_lc end
true
```

See also [`hc_triangulation`](@ref).

Role: invariant
"""
function Z_TuraevViro(hc::HeegaardCode,requests::AbstractVector;joining_disk=nothing,
        disk_name=nothing,disk_id=nothing,certificate::Bool=false,kwargs...)
    t=hc_triangulation(hc;joining_disk=joining_disk,disk_name=disk_name,disk_id=disk_id,signature=false)
    answer=Z_TuraevViro(t,requests;certificate=certificate,kwargs...)
    certificate && (answer.certificate["joining_disk"]=t.joining_disk)
    answer
end

"""
    Z_TuraevViro(inputs::AbstractVector{<:HCClosedTriangulation}, requests::AbstractVector; simplify=true, certificate=false, runner=nothing, timeout=120) -> Union{Vector,NamedTuple}

Return exact Turaev–Viro values for several finite triangulations and parameter pairs.

# Arguments
Both `inputs` and `requests` must be nonempty, checked with `ArgumentError`.
Each triangulation and each `(r, parity)` pair follows the corresponding
single-input overload. Invalid pairs raise `ArgumentError`; backend failures
propagate, and wrong input or request counts raise `ReginaTVError(:count)`.
Remaining keywords have the single-input batch meanings.

# Returns
One request-ordered value vector per input, in input order, or
`(values, certificate)` with this nested vector and the complete backend reply.

# Guarantees
1. Each result retains the backend's exact coefficients for the corresponding
   triangulation/request pair. [proved/agent: by construction]
2. Input records and request containers are not modified.
   [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
The single-input geometric and normalization limits apply. One timeout
covers the entire collection rather than each input independently.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> try Z_TuraevViro(HCClosedTriangulation[], [(3,true)]) catch e; e isa ArgumentError end
true
```

Role: invariant
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

"""
    Z_TuraevViro(inputs::AbstractVector{<:HeegaardCode}, requests::AbstractVector; certificate=false, kwargs...) -> Union{Vector,NamedTuple}

Return batches of exact Turaev–Viro values for the forward triangulations of several HCs.

# Arguments
Each HC must satisfy `hc_triangulation` with its default disk selection;
conversion errors propagate. `requests` and remaining keywords follow the
triangulation-collection overload. Empty inputs or requests and invalid
pairs raise `ArgumentError` there. This overload does not accept per-input
disk selectors; forwarded keywords are backend-evaluation keywords.

# Returns
One request-ordered vector per input, or `(values, certificate)`. When
requested, each backend job's certificate includes its selected joining disk.

# Guarantees
1. The values evaluate each input's forward triangulation in input and
   request order. [proved/agent: by construction]
2. Input HCs and the request vector are not modified.
   [proved/agent: by construction]

# Conventions
Turaev–Viro normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#moo-tv-normalization).

# Not guaranteed
Each conversion retains its own geometric evidence limits. The backend's
normalization and collection-wide timeout limitations apply.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> try Z_TuraevViro(HeegaardCode[], [(3,true)]) catch e; e isa ArgumentError end
true
```

Role: invariant
"""
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
