"""Draw checked Lickorish Data through the owner-approved placement/finish.

Input checks: research threads/PJ260922/code/ld_check.py; the .ld parser and
section 3.1 checks are reused in lickorish_input.py. K constants and face counts
are obtained from VirtualLink and lickorish_reference via the Julia worker.
No legacy marked-document layout or mesh optimizer is used.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

import lickorish_input as ld
from lickorish_placement import place
from lickorish_render import write_svg, export_png

DEBUG_MODE = False


class NativeChecker:
    """Context-managed native worker; reuse one instance for a batch of figures."""
    def __init__(self, julia='julia', project=None, compiled_modules='yes'):
        repo=Path(__file__).resolve().parents[2]
        command=[julia,'--startup-file=no',f'--compiled-modules={compiled_modules}',
                 f'--project={project or repo}',str(Path(__file__).with_name('lickorish_native_check.jl'))]
        self.process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                                      text=True,encoding='utf-8')
        ready=self.process.stdout.readline().strip()
        if ready!='READY':
            self.close()
            raise RuntimeError(f'Native checker did not start: {ready!r}')

    def __enter__(self): return self
    def __exit__(self,*_): self.close()

    def close(self):
        if self.process.stdin and not self.process.stdin.closed: self.process.stdin.close()
        self.process.wait()

    def __call__(self,data):
        self.process.stdin.write(json.dumps(data,ensure_ascii=True)+'\n');self.process.stdin.flush()
        line=self.process.stdout.readline()
        if not line: raise RuntimeError('Native checker exited without a result')
        result=json.loads(line)
        if 'error' in result: raise ValueError(result['error'])
        return result


def read_ld(path):
    path=Path(path)
    if path.suffix.lower()=='.json':
        data=json.loads(path.read_text(encoding='utf-8'))
        for key in ['places','disks','points','weights']: data.setdefault(key,{})
    else:
        g,words,places,disks,points,weights=ld.parse(path.read_text(encoding='utf-8-sig'))
        data=dict(genus=g,words=words,places=places,disks=disks,points=points,weights=weights)
    # mock3.py expects each crossing's components in input component order.
    crossings={}
    for name,word in data['words'].items():
        for label in word: crossings.setdefault(str(abs(label)),[]).append(name)
    data['crossings']=crossings
    return data


def validate(data,native):
    """Only section 4 input checks; no output quality or replay gates."""
    g,words=data['genus'],data['words']
    if not isinstance(g,int) or isinstance(g,bool) or g<0:
        raise ValueError('Genus must be a nonnegative integer')
    problems,_=ld.scan_labels(words.items())
    names=ld.k_names(g)
    if list(words)[:len(names)]!=names: problems.append(f'K components must come first: {names}')
    if problems: raise ValueError({'problems':problems,'words':words})
    report=native(data)
    for name,expected in report['constants'].items():
        actual=[x for x in words[name] if abs(x)<=3*g-2]
        if ld.canon(actual)!=ld.canon(expected): problems.append(f'{name}: {actual}; expected cyclic {expected}')
    if any(words.values()):
        if report['chi']!=2-2*g: problems.append(f"Euler characteristic {report['chi']}, expected {2-2*g}; native regions {report['regions']}")
        if report['components']!=1: problems.append(f"not connected: {report['components']} components")
    problems.extend(ld.check_extensions(g,words,data['places'],data['disks'],data['points'],data['weights']))
    if problems: raise ValueError({'problems':problems,'words':words})
    return report


def draw(data,output,native,*,steps=300,node='node',png=False,inkscape='inkscape'):
    """Validate, place, optionally finish; always preserve authoritative placement."""
    report=validate(data,native)
    output=Path(output)
    try:
        placed,_,details=place(data)
    except (AssertionError,ValueError,IndexError,StopIteration) as error:
        # Preserve the failed face/arc from mock3.py and the exact input words.
        # Placement failure is terminal; finishing is never a recovery route.
        raise ValueError({'placement_error':str(error),'words':data['words']}) from error
    placement_path=output.with_suffix('.placement.json')
    placement_path.write_text(json.dumps(placed),encoding='utf-8')
    write_svg(placed,output.with_suffix('.placement.svg'))
    result=placed
    if steps:
        target=output.with_suffix('.finished.json')
        proc=subprocess.run([node,str(Path(__file__).with_name('lickorish_relaxation.cjs')),
                            str(placement_path),str(target),str(steps)],capture_output=True,text=True)
        # Section 6.3(5): an optional finishing failure returns placement.
        if proc.returncode==0: result=json.loads(target.read_text(encoding='utf-8'))
        else: print('Finishing failed; returning placement:',proc.stderr)
    write_svg(result,output)
    if png: export_png(output,output.with_suffix('.png'),inkscape=inkscape)
    return report,details


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',nargs='+');parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--julia',default=os.environ.get('JULIA','julia'))
    parser.add_argument('--project');parser.add_argument('--compiled-modules',default='yes')
    parser.add_argument('--steps',type=int,default=300);parser.add_argument('--png',action='store_true')
    parser.add_argument('--inkscape',default='inkscape')
    args=parser.parse_args();args.output_dir.mkdir(exist_ok=True,parents=True)
    with NativeChecker(args.julia,args.project,args.compiled_modules) as native:
        for path in args.input:
            draw(read_ld(path),args.output_dir/(Path(path).stem+'.svg'),native,
                 steps=args.steps,png=args.png,inkscape=args.inkscape)


if __name__=='__main__': main()
