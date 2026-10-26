'use strict';
const $=s=>document.querySelector(s), betaColors=['#2459a9','#27806b','#8650a1','#986d15'];
let current=null,token=0,urls=[];
let currentView='alpha_cut';
function setView(view){
 if(!current)return;
 const cut=view==='alpha_cut';
 // No client-side routing or atlas cropping. Both SVGs come from Julia.
 const svg=cut?current.alpha_cut_svg:current.surface_svg;
 if(!svg)throw Error('この保存図は旧形式です。組合せデータだけで再生成してください。');
 currentView=view;$('#drawing').innerHTML=svg;$('#comparison').hidden=true;
 $('#cut').setAttribute('aria-pressed',String(cut));$('#atlas').setAttribute('aria-pressed',String(!cut));
 $('#figure-label').textContent=cut?`ALPHA-CUT · ${2*current.genus} PUNCTURE DISKS`:'参考：PUNCTURED SPHERE + ANNULI';
 download('#download-svg',svg,'image/svg+xml',current.id+(cut?'-alpha-cut.svg':'-surface.svg'));
 clearSelection();
}
function pointLabel(id){return current.marked?current.koda?.ograph.metadata.marked_projection?.intersection_labels[String(id)]||String(id):current.chart.intersections.find(p=>p.id===id)?.label||String(id);}
function download(selector,data,type,name){const u=URL.createObjectURL(new Blob([data],{type}));urls.push(u);$(selector).href=u;$(selector).download=name;}
function clearSelection(){document.querySelectorAll('.dim,.selected').forEach(e=>e.classList.remove('dim','selected'));document.querySelectorAll('#legend button,#cycles button').forEach(e=>e.setAttribute('aria-pressed','false'));$('#status').textContent=current ? `${current.intersections} 個の指定交点 · 同じラベルでα⁺・α⁻を貼り合わせ` : '';}
function selectCurve(family,id){clearSelection();document.querySelectorAll('#drawing [data-curve-family]').forEach(e=>e.classList.toggle('dim',e.dataset.curveFamily!==family||Number(e.dataset.curveId)!==id));document.querySelectorAll('#legend button').forEach(e=>e.setAttribute('aria-pressed',String(e.dataset.family===family&&Number(e.dataset.id)===id)));$('#status').textContent=`${family==='alpha'?'α':'β'}${id} を表示中`;}
function selectPoint(id){if(current.marked){const label=pointLabel(id);document.querySelectorAll('#drawing [data-crossing]').forEach(e=>e.classList.toggle('selected',e.dataset.crossing===label));$('#status').textContent=`元の図式の交点 ${label}`;return;}clearSelection();const p=current.chart.intersections.find(p=>p.id===id);if(!p)return;document.querySelectorAll('#drawing [data-link-id],#cycles [data-link-id]').forEach(e=>{const yes=Number(e.dataset.linkId)===id;e.classList.toggle('selected',yes);if(e.tagName==='BUTTON')e.setAttribute('aria-pressed',String(yes));});$('#status').textContent=`交点 ${pointLabel(id)}${p.sign===1?'+':'−'} · α${p.alpha_id} の ${p.alpha_position} 番目 · β${p.beta_id} の ${p.beta_position} 番目 · ○はこの交点へ続く接続口`;}
function showDisk(data){const d=data.chart.joining_disk;$('#joining-section').hidden=!d;$('#toggle-disk').hidden=!d;$('#toggle-disk').setAttribute('aria-pressed','true');$('#drawing').classList.remove('hide-joining-disk');$('#joining-boundary').replaceChildren();$('#joining-note').textContent=data.joining_disk_note||'';$('#joining-note').hidden=!data.joining_disk_note;if(!d)return;$('#joining-description').textContent=`領域 ${d.disk.region} · 境界 ${d.boundary.length} 辺 · ${data.visual_review==='owner-confirmed'?'見た目はオーナー確認済み':'見た目はオーナー確認待ち'}`;
for(const [k,b] of d.boundary.entries()){const button=document.createElement('button');button.textContent=pointLabel(b.from);button.onclick=()=>selectPoint(b.from);$('#joining-boundary').append(button,document.createTextNode(` → ${b.family==='alpha'?'α':'β'}${b.curve_id}${b.dart%2?'':'⁻¹'} → `));if(k===d.boundary.length-1)$('#joining-boundary').append(document.createTextNode(pointLabel(b.to)));}}
function show(data){current=data;urls.forEach(URL.revokeObjectURL);urls=[];$('#joining-section small').hidden=false;$('#figure-label').textContent='SURFACE ATLAS · ○ 接続口 · ● 数学的交点';$('#error').hidden=true;$('#title').textContent=data.title;$('#subtitle').textContent=`genus ${data.genus} · ${2*data.genus} puncture disks · ${data.chart.alpha_cycles.length} α curves · ${data.beta} β curves · ${data.intersections} signed intersections`;$('#drawing').innerHTML=data.svg;$('#base-drawing').innerHTML=data.base_svg;$('#genus').textContent=data.genus;$('#original-link').href=`/d_${data.id}.html`;
$('#legend').replaceChildren();$('#cycles').replaceChildren();
for(const family of ['alpha','beta'])for(const c of data.chart[family+'_cycles']){const label=(family==='alpha'?'α':'β')+c.id,b=document.createElement('button');b.textContent=label;b.dataset.family=family;b.dataset.id=c.id;b.style.color=family==='alpha'?'#be4c45':betaColors[(c.id-1)%4];b.onclick=()=>selectCurve(family,c.id);$('#legend').append(b);const row=document.createElement('div');row.className='cycle';const name=document.createElement('button');name.className='curve-name';name.textContent=label;name.onclick=()=>selectCurve(family,c.id);row.append(name);c.visits.forEach((id,index)=>{const p=data.chart.intersections.find(p=>p.id===id),button=document.createElement('button');button.dataset.linkId=id;button.textContent=`${pointLabel(id)}${p.sign===1?'+':'−'}`;button.onclick=()=>selectPoint(id);row.append(button,document.createTextNode(index===c.visits.length-1?'↺':'→'));});$('#cycles').append(row);}
$('#drawing').onclick=e=>{const point=e.target.closest('[data-link-id]');if(point){selectPoint(Number(point.dataset.linkId));return;}const curve=e.target.closest('[data-curve-family]');if(curve)selectCurve(curve.dataset.curveFamily,Number(curve.dataset.curveId));};
showDisk(data);$('#original-link').hidden=Boolean(data.chart.joining_disk)||data.id==='imported';
$('#convert-koda').hidden=!data.chart.joining_disk;$('#koda-section').hidden=true;$('#operation-status').textContent='';
$('#geometry').textContent='保存されたβ経路を共通描画器で実現。Dは境界dartからβ経路とα境界の区間をたどります。表示切替では座標・経路・Dを変更しません。';
$('#source').textContent=`${data.source} / ${data.seed}。数学的認証・出典の状態は入力から継承しています。`;
download('#download-svg',data.svg,'image/svg+xml',data.id+'-surface.svg');download('#download-core',JSON.stringify(data.diagram,null,2),'application/json',data.id+(data.chart.joining_disk?'.punctured.hgd.json':'.hgd.json'));download('#download-layout',JSON.stringify(data.chart,null,2),'application/json',data.id+'.surface-chart.json');download('#download-feasible',JSON.stringify(data.layout||data.chart.base_layout,null,2),'application/json',data.id+'.hgd.layout.json');document.querySelectorAll('#examples button').forEach(e=>e.classList.toggle('active',e.dataset.id===data.id));setView('alpha_cut');}
function showMarkedDisplay(d,id='imported'){
 const svg=d.svg,report=d.smooth_check;
 if(!svg||!report?.no_extra_bezier_intersections||!report.normal_port_tangents)throw Error('共通描画処理の検証済み結果が必要です');
 current={marked:true,id,title:d.title||'Heegaard diagram',diagram:d};urls.forEach(URL.revokeObjectURL);urls=[];
 $('#title').textContent=current.title;
 $('#subtitle').textContent=`genus ${d.genus} · ${d.boundaries.length} puncture disks · ${d.curves.length} cut-open curves · 平面内 ${d.crossings.length}交点`;
 $('.heading .eyebrow').textContent='HEEGAARD · MARKED DATA → 2 × g';
 $('.intro').textContent=d.reference_only?'Lickorish生成曲線。l は上下の穴の対、青が m、赤が n。同じ接続口を上下反転で貼ると閉曲線になります。穴上の点は l との交点の二つのコピーです。':'αとは独立した基準曲線で切り、等径の穴を2行×genus列に配置。赤がα、青がβ、黄色が指定されたjoining disk Dです。';
 $('.aside-note').textContent='固定された埋込みを共通処理で平滑化。穴の対は上下反転で接続。';
 $('aside details').hidden=true;
 $('#drawing').innerHTML=svg;$('#figure-label').textContent=`${d.boundaries.length} PUNCTURES · MARKED CURVES${d.disk?.length?' · JOINING DISK':''}`;
 for(const id of ['atlas','original-link','convert-koda','comparison','koda-section'])$('#'+id).hidden=true;
 $('#convert-koda').hidden=!(d.source_document&&d.disk?.length);
 $('#open-import').hidden=false;$('#regenerate').hidden=false;$('#regenerate').textContent='同じ埋込みを共通処理で描画';
 $('#cut').textContent=`2 × ${d.genus}`;$('#cut').onclick=()=>{};
 $('#operation-status').textContent=`非交差部分の最小線間隔 ${report.minimum_branch_spacing_svg?.toFixed(2)??'100以上'} · 接続口以外の穴との余白 ${report.minimum_nonport_hole_clearance.toFixed(2)} SVG単位`;
 $('#joining-section').hidden=!d.disk?.length;$('#joining-description').textContent='元の境界辺と同じ曲線を使用';
 $('#joining-boundary').textContent='';
 $('#joining-section small').textContent='Dは元の指定領域です。表示側で別の領域を選んでいません。';
 $('#toggle-disk').hidden=!d.disk?.length;$('#toggle-disk').setAttribute('aria-pressed','true');$('#toggle-disk').onclick=()=>{const el=$('#drawing [data-disk]');el.style.display=el.style.display==='none'?'':'none';$('#toggle-disk').setAttribute('aria-pressed',String(el.style.display!=='none'));};
 $('#status').textContent='同じ名前の上下境界を、保存された接続口対応に従って貼る';
 $('#legend').replaceChildren();
 const reset=()=>document.querySelectorAll('#drawing [data-curve]').forEach(e=>e.style.opacity='1');
 $('#clear-selection').onclick=reset;
 for(const c of d.curves){const b=document.createElement('button');b.textContent=c.id.replace('alpha','α').replace('beta','β');b.style.color=c.color||(c.role==='alpha'?'#ba433f':'#255dba');b.onclick=()=>document.querySelectorAll('#drawing [data-curve]').forEach(e=>e.style.opacity=e.dataset.curve===c.id?'1':'.12');$('#legend').append(b);}
 $('#drawing').onclick=()=>{};
 const details=$('.statusbar + details');details.innerHTML='<summary>共通描画処理・保存</summary><p>H.plot に明示的な埋込みと境界の貼合せ対応を渡します。交点・曲線の辺語・joining disk・2×g穴配置を保ち、曲線間と穴からの余白を検査します。表示用の切断系をαとは同一視しません。</p><p>交点・指定接続口の半径20以内は接続のための例外領域です。接続口以外で穴に接する経路は受理しません。座標のない任意入力の初期配置生成や、異なるmarking間の自動対応はこの処理の対象外です。</p><a id="marked-download" download>図式＋配置を保存</a>';
 if(d.source_figure){const f=d.source_figure,p=document.createElement('p');p.id='marked-source';p.textContent=`出典：${f.title||f.cite_key}、図${f.figure}（本文 ${f.printed_page}頁／PDF ${f.pdf_page}頁）。所蔵のOCR切り出し画像を参照。画像と転記の対応情報は保存データに含まれます。`;details.insertBefore(p,details.querySelector('p'));}
 download('#marked-download',JSON.stringify(d),'application/json','MarkedDiagram.json');
 download('#download-svg',svg,'image/svg+xml',id==='KodaGenus3Auto'?'KodaCommon.svg':id+'.svg');
 $('#error').hidden=true;history.replaceState(null,'','/surface_'+encodeURIComponent(id)+'.html');
}
async function showMarkedKoda(){const r=await fetch('/surface-generated/KodaCommon.json');if(!r.ok)throw Error('共通描画結果がありません');const d=await r.json();d.title='Koda · α / β + joining disk';showMarkedDisplay(d,'KodaGenus3Auto');}
async function select(id){if(current?.marked&&id!=='KodaGenus3Auto'){location.href='/surface_'+encodeURIComponent(id)+'.html';return;}const own=++token;try{if(id==='KodaGenus3Auto'){await showMarkedKoda();return;}const r=await fetch('/surface-generated/'+encodeURIComponent(id)+'.json');if(!r.ok)throw Error('図がありません。build_full_curves.jl の生成結果を確認してください。');const data=await r.json();if(own!==token)return;if(data.schema==='marked-cut-view/1')showMarkedDisplay(data,id);else show(data);history.replaceState(null,'','/surface_'+encodeURIComponent(id)+'.html');}catch(e){if(own===token){$('#error').textContent=e.message;$('#error').hidden=false;$('#status').textContent='表示失敗';}}}
$('#clear-selection').onclick=clearSelection;
$('#cut').onclick=()=>setView('alpha_cut');
$('#atlas').onclick=()=>setView('surface');
$('#toggle-disk').onclick=()=>{const visible=$('#toggle-disk').getAttribute('aria-pressed')!=='true';$('#toggle-disk').setAttribute('aria-pressed',String(visible));$('#drawing').classList.toggle('hide-joining-disk',!visible);};
$('#open-import').onclick=()=>$('#import-dialog').showModal();$('#close-import').onclick=()=>$('#import-dialog').close();
function showKoda(data){
 current.koda=data;const g=data.ograph,c=data.complex;$('#koda-section').hidden=false;
 $('#koda-summary').textContent=`${g.orientation.length} 頂点 · ${c.edge_incidence.length} 辺 · ${c.regions.length} 領域 · ${g.metadata.o_block_count} blocks · χ=${c.checks.euler_characteristic} · 閉性 ${data.closedness.closed?'OK':'NG'}`;
 $('#koda-drawing').innerHTML=data.svg;$('#koda-vertices').replaceChildren();
 for(const [source,v] of Object.entries(g.metadata.intersection_to_vertex).sort((a,b)=>a[1]-b[1])){const b=document.createElement('button');b.textContent=`v${v} ← ${pointLabel(Number(source))} (${g.orientation[v-1]>0?'+':'−'})`;b.onclick=()=>selectPoint(Number(source));$('#koda-vertices').append(b);}
 $('#koda-edges').replaceChildren();for(const [i,e] of data.certificate.label_map.edge_paths.entries()){const row=document.createElement('p');row.textContent=`e${i+1}: ${pointLabel(e.source[0])}/${e.source[1]} → ${pointLabel(e.target[0])}/${e.target[1]} · darts [${e.source_darts.join(', ')}]`;$('#koda-edges').append(row);}
 $('#koda-regions').textContent=c.regions.map(r=>`R${r.id} ← ${r.source_cells.join(' + ')} : [${r.boundary.join(', ')}]`).join('\n');
 $('#koda-word').textContent=JSON.stringify({gauss:g.gauss,vertex_signs:g.orientation,runs:g.metadata.runs},null,2);
 download('#download-koda',JSON.stringify(data,null,2),'application/json',current.id+'.koda.json');
}
async function diagramOperation(kind){const source=current,own=++token;$('#regenerate').disabled=true;$('#convert-koda').disabled=true;$('#operation-status').textContent=kind==='convert'?'組合せデータからKoda変換中…':source.marked?'保存された埋込みを共通処理で描画中…':'配置fixtureを使わず再生成中…';try{
 const body={view:'alpha_cut',surface_payload:true,diagram:source.diagram};
 const r=await fetch('/api/'+(kind==='convert'?'convert':'render'),{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const result=await r.json();if(!r.ok)throw Error(result.error);if(own!==token)return;
 if(kind==='convert'){showKoda(result);$('#operation-status').textContent='Koda変換完了';}else if(result.schema==='marked-cut-view/1'){showMarkedDisplay(result,source.id);}else{result.title=source.title.replace(/(?: · 自動配置)+$/,'')+' · 自動配置';show(result);history.replaceState(null,'','/surface.html');}
 }catch(e){if(own===token)$('#operation-status').textContent=e.message;}finally{$('#regenerate').disabled=false;$('#convert-koda').disabled=false;}}
$('#regenerate').onclick=()=>diagramOperation('render');
$('#convert-koda').onclick=()=>current.koda?showKoda(current.koda):diagramOperation('convert');
$('#import-form').onsubmit=async e=>{e.preventDefault();const core=$('#core-file').files[0],layout=$('#layout-file').files[0];if(!core)return;$('#render').disabled=true;$('#import-status').textContent='Juliaの共通描画経路で読込中…';try{const input=JSON.parse(await core.text()),project=input.schema==='akaghef.heegaard.project/1';const body={view:'alpha_cut',surface_payload:true,diagram:project?input.diagram:input,layout:layout?JSON.parse(await layout.text()):project?input.layout:null};const r=await fetch('/api/render',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});const data=await r.json();if(!r.ok)throw Error(data.error);++token;if(data.schema==='marked-cut-view/1')showMarkedDisplay(data);else show(data);$('#import-dialog').close();history.replaceState(null,'','/surface.html');}catch(e){$('#import-status').textContent=e.message;}finally{$('#render').disabled=false;}};
(async()=>{try{const r=await fetch('/surface-generated/manifest.json');if(!r.ok)throw Error('図の生成結果がありません。');const entries=await r.json();for(const item of entries){const b=document.createElement('button');b.dataset.id=item.id;b.textContent=item.title;b.onclick=()=>select(item.id);$('#examples').append(b);}const requested=location.pathname.match(/^\/surface_(.+)\.html$/)?.[1];const oldKoda=new URLSearchParams(location.search).get('cellulated')==='KodaGenus3';await select(requested||(oldKoda?'KodaGenus3Auto':'KodaGenus3Auto'));}catch(e){$('#error').textContent=e.message;$('#error').hidden=false;}})();
