/** cf:docs/heegaard/261109_smoothing_gates.md — existing corpus inlet and unchanged valid sampling */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const {createRelaxer}=require('../scripts/heegaard/lc_relaxation.cjs');
const [folder,beforeFile,target]=process.argv.slice(2),before=require(path.resolve(beforeFile)).createRelaxer;
const rows=JSON.parse(fs.readFileSync(path.join(folder,'results.json'))),result=[];
for(const row of rows){
 const data=JSON.parse(fs.readFileSync(path.join(folder,row.artifact.replace(/\.svg$/,'.json'))));
 const r={path:row.path};
 try{
  const model=createRelaxer(data,{DEBUG_MODE:true}),old=before(data),a=model.snapshot(),b=old.snapshot();
  r.unchanged=Object.keys(a.curves).every(k=>JSON.stringify(a.curves[k].points)===JSON.stringify(b.curves[k].points));
  r.status=model.status();assert(r.status.valid);
  r.before=old.status();r.changedCurves=[];
  for(const k of Object.keys(a.curves)){
    const p=a.curves[k].points,q=b.curves[k].points;if(JSON.stringify(p)===JSON.stringify(q))continue;
    let i=0;while(i<Math.min(p.length,q.length)&&JSON.stringify(p[i])===JSON.stringify(q[i]))i++;
    r.changedCurves.push({curve:k,beforeCount:q.length,afterCount:p.length,firstIndex:i,before:q[i]||null,after:p[i]||null});
  }
  if(!r.unchanged)assert(r.status.initialRepairs>0||r.status.inletAddedKnots>0);
  if(r.unchanged){const gate=model.debugGate(b);assert(!gate.segments.length&&!gate.vertices.length);}
 }catch(error){r.error=error.message;}
 result.push(r);fs.writeFileSync(target,JSON.stringify(result,null,2));
 console.log(row.path,r.error||JSON.stringify({unchanged:r.unchanged,repairs:r.status.initialRepairs,points:r.status.points}));
}
assert.equal(result.length,99);assert(result.every(r=>!r.error));
