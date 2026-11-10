/** cf:docs/heegaard/261109_cycle_knots.md 配置から渡す交点 — persisted knot visits survive reload and explicit resampling */
const fs=require('fs'),path=require('path');
const {createRelaxer}=require('../scripts/heegaard/lc_relaxation.cjs');
const root=process.argv[2];const receipts=[];
for(const row of JSON.parse(fs.readFileSync(path.join(root,'smoothing.json'),'utf8'))){
 const data=JSON.parse(fs.readFileSync(path.join(root,'corpus',row.artifact.replace(/\.svg$/,'.json')),'utf8'));
 if(!data.crossing_knots)continue;
 const model=createRelaxer(data),before=model.snapshot();model.resample();const after=model.snapshot();
 let error=0,sharedError=0;
 for(let i=0;i<after.crossing_knots.length;i++){
  const c=after.crossing_knots[i];const ps=c.visits.map(v=>after.curves[v.curve].points[v.index]);
  c.visits.forEach((v,j)=>{const w=before.crossing_knots[i].visits[j];const p=before.curves[w.curve].points[w.index];
   error=Math.max(error,Math.hypot(p[0]-ps[j][0],p[1]-ps[j][1]));sharedError=Math.max(sharedError,Math.hypot(ps[j][0]-ps[0][0],ps[j][1]-ps[0][1]));});
 }
 if(error>1e-7||sharedError>1e-7||model.status().crossings!==row.expected_crossings)throw Error('knot persistence failure: '+row.path);
 receipts.push({path:row.path,resamplePositionError:error,sharedKnotError:sharedError,crossings:model.status().crossings});
}
fs.writeFileSync(path.join(root,'metadata-reload.json'),JSON.stringify(receipts,null,2));console.log(JSON.stringify(receipts,null,2));
