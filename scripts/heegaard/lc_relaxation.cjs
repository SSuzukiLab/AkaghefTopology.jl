/** Numerical kernel from owner-approved app/index.html:
 * initial, resamplePiece, knotSpacing, rebuild, forces, project, stepOnce,
 * adaptiveResample. See docs/specs/lickorish-code.md section 6.3(5).
 * DOM/UI excluded; force arithmetic and defaults preserved.
 * cf:docs/heegaard/261109_cycle_knots.md 配置から渡す交点 — explicit LC knots at the inlet and persisted visit indices
 * cf:docs/heegaard/261109_smoothing_gates.md — LC crossing, local resampling, and boundary gates
 */
'use strict';
function createRelaxer(data, parameters={}) {
const par = Object.assign({"kt": 0.4, "kb": 0.1, "kr": 3.0, "rc": 220.0, "clear": 14.0, "hw": 2.0, "cw": 0.5, "rx": 40.0, "mu": 0.8, "h": 0.5, "smax": 40.0, "smin": 8.0, "spf": 3.0}, {family:false}, parameters);
let resampleRejected=false;
let DATA,W,H,holes,names,closed,raw;   // input
let X,N,S,M,off,curves,shared,fixed,coinc,same,share,portseg,base,V,posIn;
let knots;  // per curve: Set of "x,y" strings for crossing knots
let running=false,accepted=0,rejected=0,frames=0;
let crossingSigns=[],initialRepairs=0,inletAddedKnots=0,lastFailure=null;

function vsub(a,b){return [a[0]-b[0],a[1]-b[1]]}
function dist(a,b){return Math.hypot(a[0]-b[0],a[1]-b[1])}
function segInt(a,b,c,d){const r=vsub(b,a),s=vsub(d,c);const den=r[0]*s[1]-r[1]*s[0];if(Math.abs(den)<1e-12)return null;
  const q=vsub(c,a);const u=(q[0]*s[1]-q[1]*s[0])/den,v=(q[0]*r[1]-q[1]*r[0])/den;return (u>1e-9&&u<1-1e-9&&v>1e-9&&v<1-1e-9)?[u,v]:null}
function key(p){return p[0].toFixed(3)+','+p[1].toFixed(3)}

function densify(pts,cl,sp){const out=[pts[0].slice()];const n=pts.length;
  for(let i=0;i<(cl?n:n-1);i++){const a=pts[i],b=pts[(i+1)%n];const m=Math.max(1,Math.ceil(dist(a,b)/sp));
    for(let k=1;k<=m;k++){if(cl&&i===n-1&&k===m)break;out.push([a[0]+(b[0]-a[0])*k/m,a[1]+(b[1]-a[1])*k/m]);}}
  return out}

function initial(){
  // densify, insert crossing knots between different curves
  const cv={};for(const k of names)cv[k]=densify(raw[k],closed[k],12);
  const ins={};for(const k of names)ins[k]=[];
  for(let ia=0;ia<names.length;ia++)for(let ib=ia+1;ib<names.length;ib++){const a=names[ia],b=names[ib],A=cv[a],B=cv[b];
    for(let i=0;i<(closed[a]?A.length:A.length-1);i++)for(let j=0;j<(closed[b]?B.length:B.length-1);j++){
      const r=segInt(A[i],A[(i+1)%A.length],B[j],B[(j+1)%B.length]);if(r){const p=[A[i][0]+(A[(i+1)%A.length][0]-A[i][0])*r[0],A[i][1]+(A[(i+1)%A.length][1]-A[i][1])*r[0]];
        ins[a].push([i,r[0],p]);ins[b].push([j,r[1],p]);}}}
  knots={};
  for(const k of names){const Pk=cv[k].slice();ins[k].sort((x,y)=>y[0]-x[0]||y[1]-x[1]);for(const [i,u,p] of ins[k])Pk.splice(i+1,0,p.slice());cv[k]=Pk;knots[k]=new Set(ins[k].map(t=>key(t[2])));}
  // Register LC vertices already present in the polylines.
  for(const crossing of DATA.crossing_knots||[]){
    for(const visit of crossing.visits)knots[visit.curve].add(key(raw[visit.curve][visit.index]));
  }
  // A chord crossing can land exactly on a densified vertex as well. Complete
  // only visits absent from the reference strict-interior intersection pass.
  if(DATA.crossing_structure){
    const missing=k=>knots[k].size<(DATA.crossing_structure[k]||[]).length;
    const extra=Object.fromEntries(names.map(k=>[k,[]]));
    for(let ia=0;ia<names.length;ia++)for(let ib=ia+1;ib<names.length;ib++){
      const a=names[ia],b=names[ib];if(!missing(a)&&!missing(b))continue;
      const A=cv[a],B=cv[b];
      for(let i=0;i<(closed[a]?A.length:A.length-1);i++)for(let j=0;j<(closed[b]?B.length:B.length-1);j++){
        const r=vsub(A[(i+1)%A.length],A[i]),s=vsub(B[(j+1)%B.length],B[j]),q=vsub(B[j],A[i]),den=r[0]*s[1]-r[1]*s[0];
        if(Math.abs(den)<1e-12)continue;
        const u=(q[0]*s[1]-q[1]*s[0])/den,v=(q[0]*r[1]-q[1]*r[0])/den;
        if(u<-1e-10||u>1+1e-10||v<-1e-10||v>1+1e-10)continue;
        const p=[A[i][0]+r[0]*u,A[i][1]+r[1]*u];
        for(const [k,index,t] of [[a,i,u],[b,j,v]])if(!knots[k].has(key(p))){
          knots[k].add(key(p));inletAddedKnots++;
          if(t>1e-9&&t<1-1e-9)extra[k].push([index,t,p]);
        }
      }
    }
    for(const k of names){extra[k].sort((x,y)=>y[0]-x[0]||y[1]-x[1]);for(const [i,u,p] of extra[k])cv[k].splice(i+1,0,p.slice());}
  }
  return cv}

function resamplePiece(pts,sp){ // sp: target spacing per input vertex
  const L=[0];for(let i=1;i<pts.length;i++)L.push(L[i-1]+dist(pts[i],pts[i-1]));const T=L[L.length-1];
  if(T<1e-9)return [pts[0].slice(),pts[pts.length-1].slice()];
  const interp=(s,arr)=>{let i=1;while(i<L.length-1&&L[i]<s)i++;const t=(s-L[i-1])/Math.max(L[i]-L[i-1],1e-12);return arr[i-1]+(arr[i]-arr[i-1])*t};
  const cap=T/3;   // at least three segments between consecutive knots
  const out=[0];let s=0;while(true){const spv=Math.min(cap,interp(s,sp));s+=spv;if(s>=T-0.5*spv)break;out.push(s);}out.push(T);
  const result=out.map(t=>[interp(t,pts.map(p=>p[0])),interp(t,pts.map(p=>p[1]))]);
  result.stations=out;return result}

function trackedPiece(pts,sp,indices,keep){
  const cuts=[0];for(let i=1;i<pts.length-1;i++)if(keep.has(indices[i]))cuts.push(i);cuts.push(pts.length-1);
  const result=[],edges=[];
  for(let j=0;j<cuts.length-1;j++){
    const a=cuts[j],b=cuts[j+1],part=pts.slice(a,b+1),rp=resamplePiece(part,sp.slice(a,b+1));
    const L=[0];for(let i=1;i<part.length;i++)L.push(L[i-1]+dist(part[i],part[i-1]));
    const stations=rp.stations||[0,L[L.length-1]];
    for(let i=0;i<rp.length-1;i++){
      const skipped=[];for(let t=1;t<part.length-1;t++)if(L[t]>stations[i]+1e-10&&L[t]<stations[i+1]-1e-10)skipped.push(indices[a+t]);
      edges.push(skipped);
    }
    result.push(...(j?rp.slice(1):rp));
  }
  result.edges=edges;return result;
}

function knotSpacing(Pk,k,sp){
  // tighten spacing within 3*smin (arc length) of every knot and of the ends of open curves
  const n=Pk.length;const L=[0];for(let i=1;i<n;i++)L.push(L[i-1]+dist(Pk[i],Pk[i-1]));const T=L[n-1];
  const marks=[];for(let i=0;i<n;i++)if(knots[k].has(key(Pk[i])))marks.push(L[i]);if(!closed[k]){marks.push(0);marks.push(T);}
  const R=3*par.smin;const out=sp.slice();
  for(let i=0;i<n;i++)for(const m of marks){let d=Math.abs(L[i]-m);if(closed[k])d=Math.min(d,T-d);if(d<R)out[i]=Math.min(out[i],par.smin+(par.smax-par.smin)*(d/R)**2*0);}
  return out}
function rebuild(cv,spOf,retained={}){
  const nw={},provenance={};
  for(const k of names){let Pk=cv[k];const n=Pk.length;let sp=spOf?spOf[k]:new Array(n).fill(par.smax);sp=knotSpacing(Pk,k,sp);
    const original=Array.from({length:n},(_,i)=>i),keep=retained[k]||new Set();provenance[k]=[];
    let kn=[];for(let i=0;i<n;i++)if(knots[k].has(key(Pk[i])))kn.push(i);
    if(closed[k]){if(!kn.length)kn=[0];const r0=kn[0];Pk=Pk.slice(r0).concat(Pk.slice(0,r0));sp=sp.slice(r0).concat(sp.slice(0,r0));kn=kn.map(i=>(i-r0+n)%n).sort((a,b)=>a-b);kn.push(n);
      const pieces=[];for(let q=0;q<kn.length-1;q++){const a=kn[q],b=kn[q+1];const idx=[];for(let i=a;i<b;i++)idx.push(i);idx.push(b%n);const rp=trackedPiece(idx.map(i=>Pk[i]),idx.map(i=>sp[i]),idx.map(i=>(i+r0)%n),keep);provenance[k].push(...rp.edges);rp.pop();pieces.push(rp);}
      nw[k]=[].concat(...pieces);}
    else{kn=[0].concat(kn.filter(i=>i>0&&i<n-1)).concat([n-1]);const pieces=[];for(let q=0;q<kn.length-1;q++){const a=kn[q],b=kn[q+1];const rp=trackedPiece(Pk.slice(a,b+1),sp.slice(a,b+1),original.slice(a,b+1),keep);provenance[k].push(...rp.edges);pieces.push(q?rp.slice(1):rp);}nw[k]=[].concat(...pieces);}
  }
  curves=nw;off={};const pts=[];for(const k of names){off[k]=pts.length;for(const p of nw[k])pts.push(p);}
  X=pts.map(p=>[p[0],p[1]]);N=X.length;V=X.map(()=>[0,0]);
  const groups=new Map();for(let g=0;g<N;g++){const kk=key(X[g]);if(!groups.has(kk))groups.set(kk,[]);groups.get(kk).push(g);}
  shared=[...groups.values()].filter(v=>v.length>1);
  crossingSigns=[];
  if(DATA.crossing_structure){
    const labelled=new Map();
    for(const k of names){const visits=DATA.crossing_structure[k]||[];let j=0;
      nw[k].forEach((p,i)=>{if(knots[k].has(key(p))){const label=visits[j++];if(label===undefined)throw Error('LC gate: extra knot on '+k);const id=Math.abs(label);if(!labelled.has(id))labelled.set(id,[]);labelled.get(id).push({g:off[k]+i,sign:Math.sign(label)});}});
      if(j!==visits.length)throw Error('LC gate: missing knot on '+k+'; word='+JSON.stringify(visits));
    }
    const ordered=[...labelled.values()].sort((a,b)=>a[0].g-b[0].g);
    shared=ordered.map(v=>v.map(x=>x.g));crossingSigns=ordered.map(v=>v[0].sign);
  }
  S=[];for(const k of names){const n=nw[k].length,o=off[k];for(let i=0;i<(closed[k]?n:n-1);i++)S.push([o+i,o+(i+1)%n]);}M=S.length;
  fixed=new Uint8Array(N);for(const k of names)if(!closed[k]){fixed[off[k]]=1;fixed[off[k]+nw[k].length-1]=1;}
  posIn=new Array(N);for(const k of names)for(let i=0;i<nw[k].length;i++)posIn[off[k]+i]=[k,i];
  coinc=new Array(N);for(let g=0;g<N;g++)coinc[g]=[g];for(const v of shared)for(const g of v)coinc[g]=v;
  same=new Uint8Array(N*M);share=new Uint8Array(M*M);const segsOf=new Array(N);for(let g=0;g<N;g++)segsOf[g]=[];
  S.forEach(([a,b],si)=>{segsOf[a].push(si);segsOf[b].push(si);});
  S.forEach(([a,b],si)=>{for(const g of [a,b])for(const gg of coinc[g]){const [k,i]=posIn[gg];const n=nw[k].length,o=off[k];
      for(const d of [-2,-1,0,1,2]){const j=closed[k]?((i+d)%n+n)%n:i+d;if(j>=0&&j<n)same[(o+j)*M+si]=1;}
      for(const sj of segsOf[gg]){share[si*M+sj]=1;share[sj*M+si]=1;}}});
  portseg=new Uint8Array(M);S.forEach(([a,b],si)=>{if(fixed[a]||fixed[b])portseg[si]=1;});
  // All LC crossings are vertices. No additional segment crossing is allowed.
  base=new Set();
  const bad=gateBad(X),offending=new Set(bad.segments);
  for(const g of bad.vertices)S.forEach(([a,b],si)=>{if(a===g||b===g)offending.add(si);});
  if(offending.size){
    let added=0;let si=0;
    for(const k of names){if(!retained[k])retained[k]=new Set();for(const skipped of provenance[k]){
      if(offending.has(si))for(const index of skipped)if(!retained[k].has(index)){retained[k].add(index);added++;}si++;
    }}
    if(!added)throw Error('Resampling gate cannot refine: '+JSON.stringify({bad,words:DATA.crossing_structure}));
    initialRepairs+=added;rebuild(cv,spOf,retained);
  }
}
function intersections(Pt){
  const out=new Set();
  for(let i=0;i<M;i++){const a=Pt[S[i][0]],b=Pt[S[i][1]];
    for(let j=i+1;j<M;j++){const c=Pt[S[j][0]],d=Pt[S[j][1]];
      if(Math.max(a[0],b[0])+1e-9<Math.min(c[0],d[0])||Math.max(c[0],d[0])+1e-9<Math.min(a[0],b[0])||Math.max(a[1],b[1])+1e-9<Math.min(c[1],d[1])||Math.max(c[1],d[1])+1e-9<Math.min(a[1],b[1]))continue;
      const r=vsub(b,a),s=vsub(d,c),q=vsub(c,a),den=r[0]*s[1]-r[1]*s[0];
      if(Math.abs(den)>1e-12){if(share[i*M+j])continue;const u=(q[0]*s[1]-q[1]*s[0])/den,v=(q[0]*r[1]-q[1]*r[0])/den;
        if(u>=-1e-10&&u<=1+1e-10&&v>=-1e-10&&v<=1+1e-10)out.add(i*M+j);
      }else if(Math.abs(q[0]*r[1]-q[1]*r[0])<1e-10){const rr=r[0]*r[0]+r[1]*r[1];if(rr<1e-20)continue;
        const t0=(q[0]*r[0]+q[1]*r[1])/rr,t1=t0+(s[0]*r[0]+s[1]*r[1])/rr;
        if(Math.min(1,Math.max(t0,t1))-Math.max(0,Math.min(t0,t1))>1e-10)out.add(i*M+j);
      }}}
  return out}
function setEq(a,b){if(a.size!==b.size)return false;for(const x of a)if(!b.has(x))return false;return true}
function holeBad(Pt,includePorts=true){const pts=[],segs=[];
  for(let g=0;g<N;g++){if(fixed[g])continue;for(const [c,r] of holes)if(dist(Pt[g],c)<r-1e-6){pts.push(g);break;}}
  for(let si=0;si<M;si++){if(!includePorts&&portseg[si])continue;const a=Pt[S[si][0]],b=Pt[S[si][1]];for(const [c,r] of holes){if(segDist(c,a,b)<r-1e-6){segs.push(si);break;}}}
  return [pts,segs]}
function boundaryBad(Pt){const out=[];if(DATA.square)for(let g=0;g<N;g++)if(Pt[g][0]<-1e-9||Pt[g][0]>W+1e-9||Pt[g][1]<-1e-9||Pt[g][1]>H+1e-9)out.push(g);return out;}
function crossingBad(Pt){const bad=[];
  shared.forEach((v,j)=>{if(v.length!==2||dist(Pt[v[0]],Pt[v[1]])>1e-7){bad.push(...v);return;}if(!crossingSigns.length)return;
    const rays=[];for(const g of v){const [k,i]=posIn[g],n=curves[k].length,o=off[k];
      for(const d of [1,-1])rays.push(vsub(Pt[o+(i+d+n)%n],Pt[g]));}
    const angles=rays.map(p=>Math.atan2(p[1],p[0])),order=[0,1,2,3].sort((a,b)=>angles[a]-angles[b]);
    const z=order.indexOf(0),cyclic=order.slice(z).concat(order.slice(0,z));const wanted=crossingSigns[j]>0?[0,2,1,3]:[0,3,1,2];
    let invalid=cyclic.some((x,i)=>x!==wanted[i]);
    for(let i=0;i<4;i++){const a=rays[order[i]],b=rays[order[(i+1)%4]],la=Math.hypot(...a),lb=Math.hypot(...b);
      if(la<1e-10||lb<1e-10||Math.abs(a[0]*b[1]-a[1]*b[0])<1e-12*la*lb)invalid=true;}
    if(invalid)bad.push(...v);
  });return bad;}
function gateBad(Pt){const pairs=[...intersections(Pt)], [hp,hs]=holeBad(Pt),bp=boundaryBad(Pt),cp=crossingBad(Pt);
  const segments=new Set(hs);for(const x of pairs){segments.add(Math.floor(x/M));segments.add(x%M);}
  return {segments:[...segments],vertices:[...new Set([...hp,...bp,...cp])],extraCrossings:pairs.length,holeSegments:hs.length,boundaryVertices:bp.length,crossingVertices:cp.length};}
function segDist(p,a,b){const ab=vsub(b,a);const L2=ab[0]*ab[0]+ab[1]*ab[1]+1e-12;let t=((p[0]-a[0])*ab[0]+(p[1]-a[1])*ab[1])/L2;t=Math.max(0,Math.min(1,t));return Math.hypot(p[0]-a[0]-ab[0]*t,p[1]-a[1]-ab[1]*t)}
function repw(d){const c=par.clear;return d<par.rc?(c/Math.max(d,1e-9))**2/Math.max(d,1e-9):0}
function forces(Pt){
  const F=Pt.map(()=>[0,0]);let dmin=Infinity;
  for(const k of names){const o=off[k],n=curves[k].length;
    for(let i=0;i<n;i++){const g=o+i;let m,p;if(closed[k]){m=o+(i-1+n)%n;p=o+(i+1)%n;}else{if(i===0||i===n-1)continue;m=g-1;p=g+1;}
      F[g][0]+=par.kt*(Pt[m][0]+Pt[p][0]-2*Pt[g][0]);F[g][1]+=par.kt*(Pt[m][1]+Pt[p][1]-2*Pt[g][1]);
      let m2,p2;if(closed[k]){m2=o+(i-2+n)%n;p2=o+(i+2)%n;}else{if(i<2||i>n-3)continue;m2=g-2;p2=g+2;}
      for(let c=0;c<2;c++)F[g][c]-=par.kb*(Pt[m2][c]-4*Pt[m][c]+6*Pt[g][c]-4*Pt[p][c]+Pt[p2][c]);}}
  // repulsion from segments (long range), scaled by segment length
  const seglen=S.map(([a,b])=>dist(Pt[a],Pt[b]));const famOf=g=>posIn[g][0][0];
  for(let g=0;g<N;g++){if(fixed[g])continue;const p=Pt[g];
    for(let si=0;si<M;si++){if(same[g*M+si])continue;if(par.family&&famOf(S[si][0])!==famOf(g))continue;
      const a=Pt[S[si][0]],b=Pt[S[si][1]];const ab=vsub(b,a);const L2=ab[0]*ab[0]+ab[1]*ab[1]+1e-12;let t=((p[0]-a[0])*ab[0]+(p[1]-a[1])*ab[1])/L2;t=Math.max(0,Math.min(1,t));
      const D=[p[0]-a[0]-ab[0]*t,p[1]-a[1]-ab[1]*t];const d=Math.hypot(D[0],D[1]);if(d<dmin)dmin=d;
      const w=repw(d)*seglen[si]/par.smax;if(w>0){F[g][0]+=par.kr*w*D[0];F[g][1]+=par.kr*w*D[1];}}
    for(const [c,r] of holes){const D=vsub(p,c);const rr=Math.hypot(D[0],D[1]);const gap=Math.max(rr-r,1e-9);const w=repw(gap)*par.hw*gap/Math.max(rr,1e-9);F[g][0]+=par.kr*w*D[0];F[g][1]+=par.kr*w*D[1];}}
  // segment midpoints repelled by punctures too (keeps chords out of the discs)
  for(let si=0;si<M;si++){if(portseg[si])continue;const a=Pt[S[si][0]],b=Pt[S[si][1]];const mid=[(a[0]+b[0])/2,(a[1]+b[1])/2];
    for(const [c,r] of holes){const D=vsub(mid,c);const rr=Math.hypot(D[0],D[1]);const gap=Math.max(rr-r,1e-9);const w=repw(gap)*par.hw*gap/Math.max(rr,1e-9)*0.5;
      for(const g of S[si]){if(fixed[g])continue;F[g][0]+=par.kr*w*D[0];F[g][1]+=par.kr*w*D[1];}}}
  // crossings repel each other
  // crossings repel each other with a slowly decaying law: linear up to the target distance Rx
  if(shared.length>1){for(let i=0;i<shared.length;i++)for(let j=0;j<shared.length;j++){if(i===j)continue;const a=Pt[shared[i][0]],b=Pt[shared[j][0]];const D=vsub(a,b);const d=Math.hypot(D[0],D[1]);
    const w=par.cw*Math.max(0,(par.rx-d)/par.rx)/Math.max(d,1e-9);for(const g of shared[i]){F[g][0]+=par.kr*w*D[0];F[g][1]+=par.kr*w*D[1];}}}
  for(const v of shared){let fx=0,fy=0;for(const g of v){fx+=F[g][0];fy+=F[g][1];}for(const g of v){F[g][0]=fx/v.length;F[g][1]=fy/v.length;}}
  for(let g=0;g<N;g++)if(fixed[g]){F[g][0]=0;F[g][1]=0;}
  return [F,dmin]}
function project(Y){for(const v of shared){let x=0,y=0;for(const g of v){x+=Y[g][0];y+=Y[g][1];}for(const g of v){Y[g][0]=x/v.length;Y[g][1]=y/v.length;}}
  for(let g=0;g<N;g++){if(fixed[g])continue;for(const [c,r] of holes){const D=vsub(Y[g],c);const rr=Math.hypot(D[0],D[1]);const wall=r+par.clear*0.5;if(rr<wall){Y[g][0]=c[0]+D[0]/Math.max(rr,1e-9)*wall;Y[g][1]=c[1]+D[1]/Math.max(rr,1e-9)*wall;}}}
  return Y}
function ok(Y){const b=gateBad(Y);return b.segments.length===0&&b.vertices.length===0}
// Preserve the reference proposal/freeze sequence when its final candidate
// satisfies the stronger LC gate. Intermediate proposals are never published.
function proposalIntersections(Pt){
  const out=new Set(),SV=shared.map(v=>Pt[v[0]]);
  for(let i=0;i<M;i++){const a=Pt[S[i][0]],b=Pt[S[i][1]];
    for(let j=i+1;j<M;j++){if(share[i*M+j])continue;const r=segInt(a,b,Pt[S[j][0]],Pt[S[j][1]]);if(!r)continue;
      const p=[a[0]+(b[0]-a[0])*r[0],a[1]+(b[1]-a[1])*r[0]];let near=false;for(const q of SV)if(dist(p,q)<=3){near=true;break;}
      if(!near)out.add(i*M+j);}}
  return out;
}
function proposalValid(Y){if(proposalIntersections(Y).size)return false;const [p,s]=holeBad(Y,false);return !p.length&&!s.length;}
function advance(strong){
  const accept=strong?ok:proposalValid,meet=strong?intersections:proposalIntersections;
  const [F,dmin]=forces(X);let mx=0;for(const f of F)mx=Math.max(mx,Math.abs(f[0]),Math.abs(f[1]));if(mx<1e-6)return false;
  let step=Math.min(par.h,0.4*dmin/mx);
  for(let tries=0;tries<6;tries++){
    let Vn=V.map((v,g)=>fixed[g]?[0,0]:[par.mu*v[0]+step*F[g][0],par.mu*v[1]+step*F[g][1]]);
    let mv=0;for(const v of Vn)mv=Math.max(mv,Math.abs(v[0]),Math.abs(v[1]));if(mv>0.45*dmin){const s=0.45*dmin/mv;Vn=Vn.map(v=>[v[0]*s,v[1]*s]);}
    let Y=project(X.map((p,g)=>[p[0]+Vn[g][0],p[1]+Vn[g][1]]));
    if(accept(Y)){X=Y;V=Vn;accepted++;return true}
    // freeze offending vertices and retry once at this step
    const offend=new Set();const nw=meet(Y);for(const x of nw)if(!base.has(x)){offend.add(S[Math.floor(x/M)][0]);offend.add(S[Math.floor(x/M)][1]);offend.add(S[x%M][0]);offend.add(S[x%M][1]);}
    for(const x of base)if(!nw.has(x)){offend.add(S[Math.floor(x/M)][0]);offend.add(S[Math.floor(x/M)][1]);offend.add(S[x%M][0]);offend.add(S[x%M][1]);}
    const [bp,bs]=holeBad(Y,strong);for(const g of bp)offend.add(g);for(const si of bs){offend.add(S[si][0]);offend.add(S[si][1]);}
    if(strong){for(const g of boundaryBad(Y))offend.add(g);
      for(const g of crossingBad(Y)){offend.add(g);const [k,i]=posIn[g],n=curves[k].length,o=off[k];offend.add(o+(i+1)%n);offend.add(o+(i-1+n)%n);}}
    for(const g of [...offend])for(const gg of coinc[g])offend.add(gg);
    for(const g of offend){Vn[g]=[0,0];V[g]=[0,0];}
    Y=project(X.map((p,g)=>[p[0]+Vn[g][0],p[1]+Vn[g][1]]));
    if(accept(Y)){X=Y;V=Vn;accepted++;return true}
    lastFailure=gateBad(Y);step/=2;V=V.map(()=>[0,0]);rejected++;
  }
  return false}
function stepOnce(){
  const saved=[X,V.map(v=>v.slice()),accepted,rejected,lastFailure];
  const proposed=advance(false);
  if(!proposed||ok(X))return proposed;
  [X,V,accepted,rejected,lastFailure]=saved;
  return advance(true);
}
function adaptiveResample(){
  // knots = current crossing coordinates; spacing from turning angle and proximity
  for(const k of names)knots[k]=new Set();for(const v of shared)for(const g of v)knots[posIn[g][0]].add(key(X[g]));
  const sp=new Array(N).fill(par.smax);
  for(const k of names){const o=off[k],n=curves[k].length;
    for(let i=0;i<n;i++){const g=o+i;let kap=0;if(closed[k]||(i>0&&i<n-1)){const m=X[o+(i-1+n)%n],p=X[o+(i+1)%n],c=X[g];const a=vsub(p,c),b=vsub(c,m);const ang=Math.abs(Math.atan2(a[0]*b[1]-a[1]*b[0],a[0]*b[0]+a[1]*b[1]));kap=ang/Math.max(0.5*(Math.hypot(a[0],a[1])+Math.hypot(b[0],b[1])),1e-9);}
      let prox=Infinity;for(let si=0;si<M;si++){if(same[g*M+si])continue;prox=Math.min(prox,segDist(X[g],X[S[si][0]],X[S[si][1]]));}
      for(const [c,r] of holes)prox=Math.min(prox,dist(X[g],c)-r);
      sp[g]=Math.max(par.smin,Math.min(par.smax,0.35/Math.max(kap,1e-9),prox*0.8));}}
  for(const k of names){const o=off[k],n=curves[k].length;for(let it=0;it<3;it++){const v=sp.slice(o,o+n);for(let i=0;i<n;i++)sp[o+i]=Math.min(v[i],0.5*(v[(i-1+n)%n]+v[(i+1)%n]));}}
  const cv={},spk={};for(const k of names){cv[k]=X.slice(off[k],off[k]+curves[k].length).map(p=>p.slice());spk[k]=sp.slice(off[k],off[k]+curves[k].length);}
  const saved=[X,curves,off,N,S,M,shared,fixed,posIn,coinc,same,share,portseg,base,V,crossingSigns];
  try{rebuild(cv,spk);}catch(error){[X,curves,off,N,S,M,shared,fixed,posIn,coinc,same,share,portseg,base,V,crossingSigns]=saved;resampleRejected=true;}
}
function totalLength(){let L=0;for(const [a,b] of S)L+=dist(X[a],X[b]);return L}

DATA=data; W=data.W; H=data.H; holes=Object.values(data.holes);
names=Object.keys(data.curves); closed={}; raw={};
for(const k of names){closed[k]=data.curves[k].closed; let pts=data.curves[k].points;
  if(closed[k]&&dist(pts[0],pts[pts.length-1])<1e-6)pts=pts.slice(0,-1);raw[k]=pts;}
rebuild(initial());
// Explicit LC visit indices follow the retained resampling vertices.
function locateSupplied(spec){
  const locations=new Map();
  for(const k of names)curves[k].forEach((p,i)=>{const id=k+'\0'+key(p);if(!locations.has(id))locations.set(id,[]);locations.get(id).push(i);});
  return spec.map(c=>({label:c.label,visits:c.visits.map(v=>({curve:v.curve,index:locations.get(v.curve+'\0'+key(v.point)).shift()}))}));
}
let supplied=DATA.crossing_knots?locateSupplied(DATA.crossing_knots.map(c=>({label:c.label,visits:c.visits.map(v=>({curve:v.curve,point:raw[v.curve][v.index]}))}))):null;
return {
  ...(parameters.DEBUG_MODE?{debugGate(document){const pts=names.flatMap(k=>document.curves[k].points);return {...gateBad(pts),positions:posIn,meshSegments:S};}}:{}),
  step: stepOnce,
  resample(){
    const positions=supplied&&supplied.map(c=>({label:c.label,visits:c.visits.map(v=>({curve:v.curve,point:X[off[v.curve]+v.index].slice()}))}));
    adaptiveResample();if(positions)supplied=locateSupplied(positions);
  },
  snapshot(){const result=JSON.parse(JSON.stringify(data));
    for(const k of names)result.curves[k].points=X.slice(off[k],off[k]+curves[k].length).map(p=>p.slice());
    if(supplied)result.crossing_knots=JSON.parse(JSON.stringify(supplied));
    return result;},
  status(){return {points:N,segments:M,crossings:shared.length,accepted,rejected,resampleRejected,initialRepairs,inletAddedKnots,lastFailure,valid:ok(X)};}
};
}
module.exports={createRelaxer};
if(require.main===module){const fs=require('fs');
  const model=createRelaxer(JSON.parse(fs.readFileSync(process.argv[2],'utf8')));
  const steps=process.argv[4]===undefined?300:Number(process.argv[4]);
  for(let i=0;i<steps;i++)if(!model.step())break;
  fs.writeFileSync(process.argv[3],JSON.stringify(model.snapshot()));
  console.log(JSON.stringify(model.status()));
}
