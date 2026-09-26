'use strict';
const CellulatedUI={
 isProject:d=>d.schema==='akaghef.heegaard.cellulated-project/1',
 mount(){
  const panel=document.createElement('section');panel.id='cellulated-info';panel.hidden=true;
  panel.innerHTML='<p>固定した貼合せ・円板全体の細分・有向辺語を実現します。再routingはしません。無交点曲線・非連結系・入れ子も同じ形式です。</p><button id="toggle-mesh" aria-pressed="true">細分を表示</button><div class="downloads"><a id="cell-project" download>図式＋配置を保存</a><a id="cell-math" download>数学データを保存</a><a id="cell-report" download>検査reportを保存</a><a id="cell-layout" download>配置を保存</a></div><details open><summary>有向辺語・貼合せ</summary><pre id="cell-words"></pre><pre id="cell-gluings"></pre></details><details><summary>領域・cut system・保証範囲</summary><pre id="cell-regions"></pre><p>認証対象は折れ線の中心線。画素・線幅・任意のBezier平滑化は対象外です。representation hashは同位類の正規形ではありません。旧Abalone/Poincaréの経路移行とKoda出力への接続は未実装です。</p></details>';
  document.querySelector('#drawing').closest('.paper').before(panel);
  const views=document.createElement('div');views.id='cell-views';
  views.innerHTML='<button id="cell-projection">前面・背面の投影</button> <button id="cell-atlas">認証済み円板atlas</button> <button id="cell-underlay" aria-pressed="true">原画像を重ねる</button><p id="cell-view-note"></p>';
  panel.prepend(views);
  $('#cell-projection').onclick=()=>this.setView(true);$('#cell-atlas').onclick=()=>this.setView(false);
  $('#cell-underlay').onclick=()=>{const b=$('#cell-underlay'),on=b.getAttribute('aria-pressed')!=='true';b.setAttribute('aria-pressed',String(on));$('#drawing').classList.toggle('hide-source-image',!on);};
  document.querySelector('#toggle-mesh').onclick=()=>{const b=document.querySelector('#toggle-mesh'),on=b.getAttribute('aria-pressed')!=='true';b.setAttribute('aria-pressed',String(on));document.querySelector('#drawing').classList.toggle('hide-mesh',!on);};
  const nav=document.createElement('nav');nav.id='cellulated-nav';nav.setAttribute('aria-label','Marked cellulation examples');document.querySelector('aside').append(nav);
  document.querySelector('#import-dialog .hint').textContent='全体細分＋辺語の新形式は配置なしでも厳密実現します。旧marked形式は実現配置が必要です。別の経路への補完は行いません。';
 },
 async menu(){
  const r=await fetch('/cellulated-generated/manifest.json');if(!r.ok)return;
  const label=document.createElement('div');label.className='eyebrow';label.textContent='MARKED CELLULATION';$('#cellulated-nav').append(label);
  for(const item of await r.json()){const b=document.createElement('button');b.textContent=item.title;b.dataset.id=item.id;b.onclick=()=>this.select(item.id);$('#cellulated-nav').append(b);}
 },
 async select(id){const own=++token;try{const r=await fetch('/cellulated-generated/'+encodeURIComponent(id)+'.json');if(!r.ok)throw Error('Cellulated example not built');const data=await r.json();if(own!==token)return;show(data);history.replaceState(null,'','/surface.html?cellulated='+encodeURIComponent(id));}catch(e){$('#error').textContent=e.message;$('#error').hidden=false;}},
 mode(on){document.body.classList.toggle('cellulated',on);$('#cellulated-info').hidden=!on;$('#drawing').classList.remove('hide-mesh');},
 setView(projection){
  const witness=current.source_projection;this.projection=projection&&!!witness;
  $('#drawing').innerHTML=this.projection?witness.svg:current.svg;
  $('#cell-projection').setAttribute('aria-pressed',String(this.projection));$('#cell-atlas').setAttribute('aria-pressed',String(!this.projection));
  $('#cell-underlay').hidden=!this.projection;$('#toggle-mesh').hidden=this.projection;
  $('#cell-view-note').textContent=this.projection?'赤＝α、青＝β、破線＝背面。手動転写PLの同一sheet内交差と、正本の全体細分への一致を検査済み。灰色の原画像は照合用。':'全三角形・全辺を有理座標で認証したatlas。同じ接続口を貼り合わせます。';
  $('#figure-label').textContent=this.projection?'MARKED SOURCE PROJECTION · 前面／背面とD':'MARKED DISK ATLAS · 指定された辺語と貼合せ';
  download('#download-svg',this.projection?witness.svg:current.svg,'image/svg+xml',current.id+(this.projection?'.projection':'')+'.svg');this.clear();
 },
 clear(){document.querySelectorAll('#drawing .dim').forEach(e=>e.classList.remove('dim'));document.querySelectorAll('#legend button').forEach(e=>e.setAttribute('aria-pressed','false'));$('#status').textContent='保存した辺列を実現 · 中心線検査済み · 再routingなし';},
 curve(id){this.clear();document.querySelectorAll('#drawing [data-curve]').forEach(e=>e.classList.toggle('dim',e.dataset.curve!==id));document.querySelectorAll('#legend button').forEach(e=>e.setAttribute('aria-pressed',String(e.dataset.curve===id)));$('#status').textContent=id+' · 保存済み辺語を全円板で選択';},
 show(data){
  this.mode(true);$('#error').hidden=true;$('#title').textContent=data.title;const doc=data.diagram,r=data.report;
  $('#subtitle').textContent=`genus ${data.genus} · ${r.surface.boundary_components} boundaries · ${doc.curves.length} curves · ${data.intersections} signed intersections`;
  $('#drawing').innerHTML=data.svg;$('#base-drawing').replaceChildren();$('#legend').replaceChildren();$('#cycles').replaceChildren();
  $('#original-link').hidden=true;$('#comparison').hidden=true;$('#convert-koda').hidden=!(r.joining_disks.length===1&&r.joining_disks[0].koda_joining);$('#koda-section').hidden=true;$('#operation-status').textContent='';
  $('#regenerate').textContent='保存した図式を再実現';$('#toggle-mesh').setAttribute('aria-pressed','true');
  $('#marking-status').textContent=`Marked whole diagram · ${doc.surface.marking_id} · representation ${r.representation_sha256}`;
  $('#mark-current').hidden=true;$('#toggle-seams').hidden=true;$('#download-marked-project').hidden=true;$('#marking-details').hidden=true;
  for(const c of doc.curves){const b=document.createElement('button');b.textContent=c.id+' ('+c.role+')';b.dataset.curve=c.id;b.onclick=()=>this.curve(c.id);$('#legend').append(b);}
  $('#drawing').onclick=e=>{const c=e.target.closest('[data-curve]');if(c)this.curve(c.dataset.curve);};
  $('#cell-words').textContent=doc.curves.map(c=>c.id+' ['+(c.closed?'closed':'arc')+']: '+c.word).join('\n');
  $('#cell-gluings').textContent=doc.surface.gluings.map(([a,b])=>a.join('/')+' ↔ '+b.join('/')+' (t ↔ 1−t)').join('\n');
  $('#cell-regions').textContent=JSON.stringify({alpha_cut:r.alpha_cut,beta_cut:r.beta_cut,regions:r.regions.map(({triangles,boundary_cycles,...rest})=>rest),joining_disks:r.joining_disks},null,2);
  $('#cell-regions').nextElementSibling.textContent='atlasの認証対象は折れ線中心線。Koda変換はcellular・complete cut systems・有効な向き付きDを満たすαβ閉曲線の入力に対応。非cellular・aux付き・複数Dの自動選択は行いません。旧図式の経路移行は別課題です。';
  const disks=r.joining_disks;$('#joining-section').hidden=!disks.length;$('#toggle-disk').hidden=!disks.length;$('#toggle-disk').setAttribute('aria-pressed','true');$('#drawing').classList.remove('hide-joining-disk');
  $('#joining-description').textContent=disks.map(d=>`${d.id}: embedded closed disk = ${d.embedded_closed_disk}; Koda joining = ${d.koda_joining}`).join(' / ');
  $('#joining-boundary').textContent=disks.map(d=>`seed ${d.seed.join(':')} · all arcs connected: ${d.arcs_all_nonempty_connected} · outward-positive: ${d.boundary_outward_positive}`).join('\n');
  $('#joining-note').hidden=false;$('#joining-note').textContent='幾何的な円板とKodaの向き条件を区別します。曲線の向きは自動反転しません。';
  $('#joining-section small').hidden=true;
  $('#figure-label').textContent='MARKED DISK ATLAS · 指定された辺語と貼合せ';
  const json=x=>JSON.stringify(x,null,2);
  download('#download-svg',data.svg,'image/svg+xml',data.id+'.svg');
  download('#cell-project',json({schema:'akaghef.heegaard.cellulated-project/1',diagram:doc,layout:data.layout}),'application/json',data.id+'.cellulated-project.json');
  download('#cell-math',json(doc),'application/json',data.id+'.cellulated.json');
  download('#cell-report',json(r),'application/json',data.id+'.report.json');
  download('#cell-layout',json(data.layout),'application/json',data.id+'.cellulated-layout.json');
  document.querySelectorAll('#examples button').forEach(b=>b.classList.remove('active'));
  document.querySelectorAll('#cellulated-nav button').forEach(b=>b.classList.toggle('active',b.dataset.id===data.id));
  $('#cell-projection').hidden=!data.source_projection;$('#cell-underlay').setAttribute('aria-pressed','true');$('#drawing').classList.remove('hide-source-image');
  $('#cell-words').closest('details').open=false;this.setView(!!data.source_projection);
 }
};
