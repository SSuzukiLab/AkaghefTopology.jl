"""Direct-driver regression: VirtualLink regions vs the independent LC checker.

Source: research threads/PJ260922/code/lc_check.py (section 4 oracle),
examples/VLExample.jl and examples/OgraphPlotSuite.jl (existing planar cases).
This is deliberately a test, never a drawing-time double implementation.
Run with research root, scratch output directory, Julia executable as arguments.
The bare-module driver does NOT replace the package test suite.
"""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def run(research, output, julia):
    research, output = Path(research), Path(output)
    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('oracle', research/'threads/PJ260922/code/lc_check.py')
    oracle = importlib.util.module_from_spec(spec); spec.loader.exec_module(oracle)
    records, checks = [], []
    for folder in ['261107-lickorish-data', '261107-ha9-curve-figures']:
        for path in sorted((research/'threads/PJ260922/data'/folder).rglob('*.lc')):
            record = dict(path=path.relative_to(research).as_posix())
            try:
                g, words, places, disks, points, weights = oracle.parse(path.read_text(encoding='utf-8-sig'))
                problems = oracle.check(g, words, places, disks, points, weights)
                record.update(genus=g, words=words, places=places, problems=problems,
                              input_status='fail' if problems else 'pass')
                if not problems:
                    active = {k:v for k,v in words.items() if v}
                    labels = sorted({abs(x) for w in active.values() for x in w})
                    compact = {x:i+1 for i,x in enumerate(labels)}
                    code = [[(1 if x>0 else -1)*compact[abs(x)] for x in w] for w in active.values()]
                    faces = oracle.surface(active)[2]
                    record.update(vertices=len(labels), faces=faces)
                    literal = 'Vector{Int}[' + ','.join('Int'+json.dumps(w) for w in code) + ']'
                    checks.append(f'check_case({json.dumps(record["path"])}, {literal}, {len(labels)}, {faces})')
            except Exception as error:
                record.update(input_status='unsupported', problems=[str(error)])
            records.append(record)
    output.mkdir(parents=True, exist_ok=True)
    (output/'261109_corpus-inputs.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    source = r'''
module Current
end
module Previous
end
function load_kernel(mod, path)
    src=read(path,String)
    stop=first(findfirst("function hmove_matrix",src))-1
    Base.include_string(mod,src[1:stop],path)
    begin_weight=first(findfirst("function set_weight!",src))
    stop_weight=first(findnext("\nend",src,begin_weight))+3
    Base.include_string(mod,src[begin_weight:stop_weight],path)
end
load_kernel(Current,joinpath(REPO,"src","virtual_link.jl"))
load_kernel(Previous,joinpath(OUT,"261109_virtual_link_before_regions.jl"))
for example in ("VLExample.jl","OgraphPlotSuite.jl")
    src=replace(read(joinpath(REPO,"examples",example),String),"using Topology"=>"")
    Base.include_string(Current,src,joinpath(REPO,"examples",example))
end
fixtures=Current.run_vl_plot_calls()
for cases in values(Current.run_ograph_plot_suite()); append!(fixtures,cases); end
for case in fixtures
    pd=Current.pd_code(case.link)
    @assert Current.regions(pd)==Previous.regions(pd) case.name
end
@assert Current.regions([1 2 1 2])==[[2,1,-2,-1]]
println("existing planar/o-graph fixtures unchanged: ",length(fixtures))
function check_case(name,code,n,expected)
    v=Current.VirtualLink()
    Current.set_data!(v;gauss=code,orientation=ones(Int,n))
    actual=Current.regions(Current.pd_code(v))
    @assert length(actual)==expected (name,expected,actual)
    @assert sum(length,actual;init=0)==4*n (name,actual)
    println("PASS\t",name,"\t",length(actual))
end
'''
    source = f'const REPO={json.dumps(repo.as_posix())}\nconst OUT={json.dumps(output.as_posix())}\n'+source+'\n'+'\n'.join(checks)
    driver=output/'261109_regions_direct_driver.jl'; driver.write_text(source, encoding='utf-8')
    proc=subprocess.run([julia,'--startup-file=no',str(driver)],capture_output=True,text=True,encoding='utf-8')
    (output/'261109_regions_direct_driver.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
    print(proc.stdout,proc.stderr)
    print('Inputs:',len(records),'eligible:',len(checks))
    return proc.returncode


if __name__=='__main__': sys.exit(run(*sys.argv[1:]))
