
const select=document.getElementById('choose');
function refresh(){const v=select.value;document.getElementById('original').style.opacity=v==='none'?'1':'.18';document.getElementById('labels').style.display=document.getElementById('showlabels').checked?'':'none';for(const p of document.querySelectorAll('.trace')){const show=(v==='all'||v===p.dataset.curve)&&(p.dataset.sheet==='F'?document.getElementById('front').checked:document.getElementById('back').checked);p.style.display=show?'':'none';}}
for(const x of document.querySelectorAll('.controls select,.controls input'))x.addEventListener('change',refresh);
for(const x of document.querySelectorAll('.point'))x.addEventListener('click',()=>{for(const r of document.querySelectorAll('.highlight'))r.classList.remove('highlight');const row=document.getElementById('row-'+x.dataset.point);row.classList.add('highlight');row.scrollIntoView({behavior:'smooth',block:'center'});});
refresh();
