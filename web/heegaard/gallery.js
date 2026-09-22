'use strict';
const $=s=>document.querySelector(s);
const names={'Source fixture':'参考図','Lens space':'レンズ空間','Stabilization':'S³・安定化','Connected sum':'連結和'};
let catalog=[],selection={representatives:[],comparisons:[]};
function homology(h){if(!h)return 'H₁ —';const terms=[];if(h.free_rank)terms.push(h.free_rank===1?'ℤ':`ℤ^${h.free_rank}`);terms.push(...h.torsion.map(n=>`ℤ/${n}`));return `H₁ = ${terms.length?terms.join(' ⊕ '):'0'}`;}
function card(e,description){
 const card=document.createElement('a');card.className='diagram-card';card.href=`/d_${encodeURIComponent(e.id)}.html`;
 const picture=document.createElement('div');picture.className='card-image';const img=document.createElement('img');img.src=`/generated/${encodeURIComponent(e.id)}.svg`;img.alt=`${e.title} のα-cut図`;img.loading='lazy';img.decoding='async';picture.append(img);
 const copy=document.createElement('div');copy.className='card-copy';const category=document.createElement('span');category.className='card-category';category.textContent=names[e.category]||'DIAGRAM';
 const title=document.createElement('h2');title.textContent=e.title;
 const meta=document.createElement('div');meta.className='card-meta';for(const value of [`genus ${e.genus}`,`${e.intersections} 交点`,`${2*e.genus} disks`]){const span=document.createElement('span');span.textContent=value;meta.append(span);}
 const note=document.createElement('p');note.className='card-description';note.textContent=description;
 const h=document.createElement('div');h.className='card-homology';h.textContent=homology(e.h1);const arrow=document.createElement('span');arrow.className='card-open';arrow.textContent='↗';h.append(arrow);
 copy.append(category,title,meta,note,h);card.append(picture,copy);return card;
}
function render(){
 const query=$('#search').value.trim().toLowerCase().normalize('NFD').replace(/[\u0300-\u036f]/g,'').replace(/\s/g,'');
 const entries=selection.representatives.map(s=>({...catalog.find(e=>e.id===s.id),description:s.description})).filter(e=>e.id);
 const items=entries.filter(e=>(!$('#genus').value||e.genus===Number($('#genus').value))&&(!query||`${e.title} ${e.id} ${e.description}`.toLowerCase().normalize('NFD').replace(/[\u0300-\u036f]/g,'').replace(/\s/g,'').includes(query)));
 $('#gallery-grid').replaceChildren(...items.map(e=>card(e,e.description)));
 $('#result-count').textContent=`代表図 ${items.length} 件`;$('#empty').hidden=items.length!==0;
 const params=new URLSearchParams();for(const id of ['search','genus'])if($('#'+id).value)params.set(id,$('#'+id).value);
 history.replaceState(null,'',location.pathname+(params.size?'?'+params.toString():''));
}
for(const id of ['search','genus'])$('#'+id).addEventListener(id==='search'?'input':'change',render);
$('#clear-filters').onclick=()=>{for(const id of ['search','genus'])$('#'+id).value='';render();};
(async()=>{try{
 const responses=await Promise.all([fetch('/generated/manifest.json'),fetch('/gallery-selection.json')]);
 if(responses.some(r=>!r.ok))throw Error('一覧を取得できませんでした。');
 [catalog,selection]=await Promise.all(responses.map(r=>r.json()));
 $('#coverage-note').textContent=selection.missing||'';$('#coverage-note').hidden=!selection.missing;
 for(const g of [...new Set(selection.representatives.map(s=>catalog.find(e=>e.id===s.id)?.genus).filter(Boolean))].sort((a,b)=>a-b)){const option=document.createElement('option');option.value=g;option.textContent=`genus ${g}`;$('#genus').append(option);}
 for(const s of selection.comparisons){const e=catalog.find(e=>e.id===s.id);if(e)$('#comparison-grid').append(card(e,s.description));}
 const selected=new Set([...selection.representatives,...selection.comparisons].map(s=>s.id));
 const groups=new Map();for(const e of catalog.filter(e=>!selected.has(e.id))){const key=`${names[e.category]||e.category} · genus ${e.genus}`;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(e);}
 for(const [title,entries] of groups){const details=document.createElement('details');const summary=document.createElement('summary');summary.textContent=title;const links=document.createElement('div');links.className='variant-links';for(const e of entries){const a=document.createElement('a');a.href=`/d_${encodeURIComponent(e.id)}.html`;a.textContent=e.title;links.append(a);}details.append(summary,links);$('#variant-groups').append(details);}
 $('#variant-groups').parentElement.hidden=groups.size===0;
 const params=new URLSearchParams(location.search);for(const id of ['search','genus'])if(params.has(id))$('#'+id).value=params.get(id);
 render();
 }catch(e){$('#gallery-error').textContent=e.message;$('#gallery-error').hidden=false;$('#result-count').textContent='読み込み失敗';}})();
