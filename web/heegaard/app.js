'use strict';
const $=s=>document.querySelector(s), colors=['#2459a9','#b75029','#27806b','#8650a1','#986d15'];
let current,manifest=[],request=0,comparison=false,selectedBeta=null,views=[],urls=[];
function download(id,data,type,name){const u=URL.createObjectURL(new Blob([data],{type}));urls.push(u);$(id).href=u;$(id).download=name;}
function error(message){$('#error').textContent=message;$('#error').hidden=!message;}
function unionViewBox(a,b){let p=a.split(/\s+/).map(Number),q=b.split(/\s+/).map(Number);let x=Math.min(p[0],q[0]),y=Math.min(p[1],q[1]);return [x,y,Math.max(p[0]+p[2],q[0]+q[2])-x,Math.max(p[1]+p[3],q[1]+q[3])-y];}
function setViews(box){views.forEach(v=>v.setAttribute('viewBox',box.join(' ')));}
function fit(){if(!views.length)return;setViews(comparison?unionViewBox(current.beforeBox,current.afterBox):current.afterBox.split(/\s+/).map(Number));}
function zoom(factor){if(!views.length)return;const b=views.at(-1).viewBox.baseVal;setViews([b.x+b.width*(1-factor)/2,b.y+b.height*(1-factor)/2,b.width*factor,b.height*factor]);}
function highlightBeta(id){selectedBeta=selectedBeta===id?null:id;document.querySelectorAll('.drawing [data-curve-id]').forEach(p=>p.style.opacity=selectedBeta===null||Number(p.dataset.curveId)===selectedBeta?'1':'.13');document.querySelectorAll('#legend button').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.beta)===selectedBeta)));}
function activateCanvas(container){const svg=container.querySelector('svg');if(!svg)return;let start=null;
container.onwheel=e=>{e.preventDefault();zoom(e.deltaY<0?.87:1.15);};
container.onpointerdown=e=>{if(e.target.hasAttribute('data-intersection-id')||e.target.hasAttribute('data-intersection-label'))return;start={x:e.clientX,y:e.clientY,box:[svg.viewBox.baseVal.x,svg.viewBox.baseVal.y,svg.viewBox.baseVal.width,svg.viewBox.baseVal.height]};container.setPointerCapture(e.pointerId);};
container.onpointermove=e=>{if(!start)return;const rect=svg.getBoundingClientRect(),scale=Math.max(start.box[2]/rect.width,start.box[3]/rect.height);setViews([start.box[0]-(e.clientX-start.x)*scale,start.box[1]-(e.clientY-start.y)*scale,start.box[2],start.box[3]]);};
container.onpointerup=()=>start=null;container.onpointercancel=()=>start=null;
container.onclick=e=>{const id=e.target.getAttribute('data-intersection-id')||e.target.getAttribute('data-intersection-label');if(!id)return;document.querySelectorAll('.drawing [data-intersection-id],.drawing [data-intersection-label]').forEach(p=>p.classList.toggle('matched',(p.dataset.intersectionId||p.dataset.intersectionLabel)===id));$('#status').textContent=`交点 ${id} の貼り合わせ相手を表示中`;};
}
function show(data){current=data;selectedBeta=null;error('');urls.forEach(URL.revokeObjectURL);urls=[];
$('#title').textContent=data.title;$('#subtitle').textContent=`genus ${data.genus}  ·  ${data.intersections} intersections  ·  ${data.layout.routes.length} β arcs  ·  ${2*data.genus} punctures`;
$('#after-svg').innerHTML=data.svg;$('#before-svg').innerHTML=data.before_svg;views=[$('#before-svg svg'),$('#after-svg svg')].filter(Boolean);data.beforeBox=views[0].getAttribute('viewBox');data.afterBox=views[1].getAttribute('viewBox');fit();
['#before-svg','#after-svg'].forEach(x=>activateCanvas($(x)));$('#legend').replaceChildren();
for(const id of [...new Set(data.layout.routes.map(r=>r.curve_id))]){const b=document.createElement('button');b.dataset.beta=id;b.setAttribute('aria-pressed','false');b.style.setProperty('--curve',colors[(id-1)%colors.length]);b.innerHTML=`<span class="swatch"></span>β${id}`;b.onclick=()=>highlightBeta(id);$('#legend').append(b);}
const m=data.layout.backend_trace.smoothing||{};$('#status').textContent=data.layout.report.valid?'接続・端点・非交差を確認済み':'未検証';
const distance=m.minimum_flattened_separation==null?'対象なし（β弧1本）':Number(m.minimum_flattened_separation).toFixed(4);
$('#metrics').textContent=`曲線間の最小距離 ${distance} · 出口は境界の法線方向（許容角15°） · ${m.accepted_steps||0} updates · ${m.stop_reason||'saved layout'}。距離は平坦化した線分での測定値。`;
$('#source').textContent=`${data.source} / ${data.seed}`;
download('#download-svg',data.svg,'image/svg+xml',`${data.id}.svg`);download('#download-core',JSON.stringify(data.diagram,null,2),'application/json',`${data.id}.hgd.json`);download('#download-layout',JSON.stringify(data.layout,null,2),'application/json',`${data.id}.hgd.layout.json`);download('#download-before',JSON.stringify(data.before,null,2),'application/json',`${data.id}-before.hgd.layout.json`);
document.querySelectorAll('#examples button').forEach(b=>b.classList.toggle('active',b.dataset.id===data.id));
}
async function select(id){const token=++request;try{$('#status').textContent='読み込み中…';const r=await fetch(`/generated/${encodeURIComponent(id)}.json`);if(!r.ok)throw Error('描画データがありません。生成処理の完了後に再読込してください。');const data=await r.json();if(token!==request)return;show(data);history.replaceState(null,'',`/d_${id}.html`);}catch(e){if(token===request)error(e.message);}}
$('#single').onclick=()=>{comparison=false;$('#canvases').classList.remove('compare');$('#single').setAttribute('aria-pressed','true');$('#compare').setAttribute('aria-pressed','false');fit();};
$('#compare').onclick=()=>{comparison=true;$('#canvases').classList.add('compare');$('#single').setAttribute('aria-pressed','false');$('#compare').setAttribute('aria-pressed','true');fit();};
$('#zoom-in').onclick=()=>zoom(.8);$('#zoom-out').onclick=()=>zoom(1.25);$('#fit').onclick=fit;
$('#open-import').onclick=()=>$('#import-dialog').showModal();$('#close-import').onclick=()=>$('#import-dialog').close();
$('#import-form').onsubmit=async e=>{e.preventDefault();const core=$('#core-file').files[0],layout=$('#layout-file').files[0];if(!core)return;$('#render').disabled=true;$('#import-status').textContent='Juliaで描画中。初回は起動・コンパイルを含みます…';try{const body={diagram:JSON.parse(await core.text()),layout:layout?JSON.parse(await layout.text()):null};const r=await fetch('/api/render',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const data=await r.json();if(!r.ok)throw Error(data.error);++request;show(data);$('#import-dialog').close();history.replaceState(null,'','/');}catch(e){$('#import-status').textContent=e.message;}finally{$('#render').disabled=false;}};
(async()=>{try{const response=await fetch('/generated/manifest.json');if(!response.ok)throw Error('図を生成しています。少し待って再読込してください。');manifest=await response.json();for(const item of manifest){const b=document.createElement('button');b.dataset.id=item.id;b.textContent=item.title;const small=document.createElement('small');small.textContent=`genus ${item.genus} · ${item.intersections} intersections`;b.append(small);b.onclick=()=>select(item.id);$('#examples').append(b);}const id=location.pathname.match(/^\/d_(.+)\.html$/)?.[1];await select(manifest.some(x=>x.id===id)?id:manifest[0].id);}catch(e){error(e.message);}})();
