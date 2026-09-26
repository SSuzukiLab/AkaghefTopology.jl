"""Build a self-contained HTML review of source trace, code, and comparison."""
import base64, json, math, html
from pathlib import Path
import build_case as B
import trace_input as T
root=B.ROOT
core=json.loads((root/'data/figure6.core.json').read_text())
koda=json.loads((root/'data/figure6.koda.json').read_text())
comparison=json.loads((root/'artifacts/comparison.json').read_text())
validation=json.loads((root/'artifacts/validation.json').read_text())
img=base64.b64encode((root/'data/source.png').read_bytes()).decode()
labelpos={'c1':(382,551),'c2':(561,408),'c3':(750,509),'c4':(739,592),'c5':(615,662),'c6':(450,700),
          'x1':(360,428),'x2':(217,564),'x3':(517,369),'x4':(887,414),'x5':(1067,474)}
svg=[f'<svg viewBox="0 0 1700 840" aria-label="原画像と独立に転写した交点・経路"><image id="original" href="data:image/png;base64,{img}" width="1700" height="840"/>',
     '<defs><marker id="dir" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="4" markerHeight="4" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10z" fill="context-stroke"/></marker></defs>']
for cid,pieces in T.PIECES.items():
    for sheet,vs in pieces:
        pts=' '.join(f'{p[0]},{p[1]}' for p in map(B.point,vs))
        role='alpha' if cid.startswith('alpha') else 'beta'
        svg.append(f'<polyline class="trace {role}" data-curve="{cid}" data-sheet="{sheet}" points="{pts}" marker-mid="url(#dir)" style="display:none"/>')
svg.append('<g id="labels">')
for name,(lx,ly) in labelpos.items():
    x,y=T.POINTS[name];text=name+(' (back)' if name=='x2' else '')
    w=132 if name=='x2' else 44
    svg.append(f'<g class="point" data-point="{name}"><path d="M{x} {y} L{lx+w/2} {ly+14}" stroke="black" stroke-width="2"/><rect x="{lx-4}" y="{ly-8}" width="{w+8}" height="35" rx="4" fill="white" stroke="black"/><text x="{lx}" y="{ly+18}" font-size="26" font-family="sans-serif">{text}</text><circle cx="{x}" cy="{y}" r="5" fill="white" stroke="black" stroke-width="2"/></g>')
svg.append('</g></svg>')
cyclesrows=''.join('<tr><th>'+html.escape(c['id'])+'</th><td><code>'+' → '.join(c['cycle'])+' → …</code></td></tr>' for family in ('alpha','beta') for c in core[family])
xsrows=''.join(f'<tr id="row-{x["id"]}"><th>{x["id"]}</th><td>{x["alpha"]}</td><td>{x["beta"]}</td><td>{x["sign"]:+d}</td><td>{"背面" if x["id"]=="x2" else "前面"}</td><td>{"Dの角 → 除去" if x["id"].startswith("c") else "o-graph頂点"}</td></tr>' for x in sorted(core['intersections'],key=lambda x:x['id']))
blocksrows=''.join(f'<tr><th>{b["curve"]}</th><td><code>{" → ".join(b["path"])}</code></td><td><code>{" ".join(("−" if b["curve"].startswith("alpha") else "+")+v[1:] for v in b["path"] if v.startswith("x"))}</code></td></tr>' for b in koda['blocks'])
def wordhtml(w):return ' '.join(f'<span class="token">{x:+d}</span>' for x in w)
actual=koda['o_data'][0][0];target=comparison['target'][0][0]
comparecells=''.join(f'<tr><th>{i+1}</th><td>{a:+d}</td><td class="{("diff" if a!=b else "")}">{b:+d}</td></tr>' for i,(a,b) in enumerate(zip(actual,target)))
# A Gauss representation of the o-graph, not a claim of a planar embedded curve.
coords=[];g=['<svg viewBox="0 0 600 550" aria-label="o-graphのGauss表示"><circle cx="300" cy="267" r="205" fill="none" stroke="#333" stroke-width="2"/>']
for j,v in enumerate(actual):
    t=-math.pi/2+2*math.pi*j/len(actual);coords.append((300+205*math.cos(t),267+205*math.sin(t)))
for n in range(1,6):
    i,j=actual.index(-n),actual.index(n);x,y=coords[i];u,v=coords[j]
    g.append(f'<line x1="{x}" y1="{y}" x2="{u}" y2="{v}" stroke="#a0a0a0" stroke-width="2"/><text x="{(x+u)/2}" y="{(y+v)/2}" font-size="19" fill="#111" stroke="white" stroke-width="5" paint-order="stroke">v{n}: {koda["vertices"][n-1]["type"]}</text>')
for j,(v,(x,y)) in enumerate(zip(actual,coords)):
    t=-math.pi/2+2*math.pi*j/len(actual)
    g.append(f'<circle cx="{x}" cy="{y}" r="4" fill="#222"/><text x="{300+231*math.cos(t)}" y="{267+231*math.sin(t)+7}" text-anchor="middle" font-size="24">{v:+d}</text>')
g.append('</svg>')
page='''<!doctype html><html lang="ja"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>種数3 Heegaard図式 — 画像転写とKoda照合</title><style>
:root{color-scheme:light}*{box-sizing:border-box}body{margin:0;background:#f5f6f7;color:#182028;font-family:system-ui,"Noto Sans CJK JP",sans-serif;line-height:1.7}main{max-width:1440px;margin:auto;padding:36px 30px}h1{font-size:30px;margin:0 0 8px}h2{font-size:23px;margin:0 0 18px}h3{font-size:18px}p{margin:10px 0 15px}.lead{font-size:18px}.badge{font-size:12px;letter-spacing:1px;color:#65717d}.card{background:white;border:1px solid #dce0e4;border-radius:10px;padding:25px;margin:22px 0}.verdict{border-left:6px solid #9b3329}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:20px}.metric{background:#eef1f3;padding:12px 16px;border-radius:6px}.metric strong{font-size:26px;display:block}.metric span{font-size:13px}.two{display:grid;grid-template-columns:1fr 1fr;gap:25px}.figure{overflow:auto}.figure>svg{width:100%;min-width:750px;height:auto}.trace{fill:none;stroke-width:6;stroke-linejoin:round;stroke-linecap:round;pointer-events:none}.trace.alpha{stroke:#be3b33}.trace.beta{stroke:#256da3}.trace[data-sheet="B"]{stroke-dasharray:13 9}.point{cursor:pointer}.controls{display:flex;flex-wrap:wrap;gap:18px;align-items:center;margin-bottom:15px}select,button{font:inherit;padding:5px 10px;border:1px solid #bbb;border-radius:5px;background:white}.note{font-size:13px;color:#52616a}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e0e4e7}thead th{background:#edf0f2}code,pre{font-family:ui-monospace,Consolas,monospace}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#eef1f3;padding:15px;line-height:1.7;border-radius:6px}.token{display:inline-block;min-width:2.6em;text-align:center;background:#eef1f3;border:1px solid #d7dde2;border-radius:4px;padding:7px 5px;margin:3px}.diff{background:#fcece8;font-weight:700}.highlight{outline:2px solid #ac3b2e;background:#fff4ed}a{color:#185d8e}footer{font-size:13px;color:#56616d;padding:20px 0}details{margin-top:20px}summary{cursor:pointer;font-weight:600}@media(max-width:800px){main{padding:18px 14px}.two{grid-template-columns:1fr}.metrics{grid-template-columns:1fr 1fr}.card{padding:16px}h1{font-size:24px}}
</style><main><div class="badge">SOURCE IMAGE → MARKED SURFACE → KODA CONVERSION</div><h1>種数3・joining disk付きHeegaard図式</h1><p class="lead">画像の前面／背面を分けて転写。提示されたo-dataは再構成の入力に使わず、最後に照合した。</p>
<section class="card verdict"><h2>結果：図式の転写は検証済み。提示o-dataとは不一致。</h2><p>画像から得た語は、Koda論文Fig. 6のE-dataと一致する。提示値は頂点4・5の両通過符号が逆であり、巡回移動・頂点番号の付け替えでは一致しない。</p><pre>computed: [-1,-2, 3, 1, 4,-5,-4, 5,-3, 2]
 supplied: [-1,-2, 3, 1,-4, 5, 4,-5,-3, 2]
    types: [-1, 1, 1, 1,-1]  (両方同じ)</pre><div class="metrics"><div class="metric"><strong>3</strong><span>曲面の種数</span></div><div class="metric"><strong>11</strong><span>真のα–β交点</span></div><div class="metric"><strong>6 → 除去</strong><span>Dの角</span></div><div class="metric"><strong>5</strong><span>o-graph頂点</span></div></div></section>
<section class="card"><h2>原画像・交点・経路の監査</h2><div class="controls"><label>転写経路 <select id="choose"><option value="none">表示しない（原画像）</option value="all">全曲線</option>'''+''.join(f'<option value="{c}">{c}</option>' for c in T.PIECES)+'''</select></label><label><input id="front" type="checkbox" checked> 前面</label><label><input id="back" type="checkbox" checked> 背面</label><label><input id="showlabels" type="checkbox" checked> 交点ラベル</label></div><div class="figure">'''+''.join(svg)+'''</div><p class="note">c1…c6はDの角。x1…x5は変換後に残る頂点。x2は背面での真の交点。前面の曲線と背面の曲線の画面上の交差は、曲面上の交点には数えない。経路を選ぶと原画像を薄くして転写折れ線を重ねる。前面／背面の切替は転写線にのみ適用する。交点ラベルを押すと下表を強調する。</p></section>
<section class="card"><h2>数学データ</h2><p>向きはすべてDを右に見る向き。交点符号εは順序付き接ベクトル(α̇,β̇)と曲面の向きを比較する。表側の画素座標ではdx∧d(−y)、裏側では逆向きで計算した。</p><div class="two"><div><h3>有向巡回順序</h3><table><thead><tr><th>曲線</th><th>交点列</th></tr></thead><tbody>'''+cyclesrows+'''</tbody></table><h3>Dの境界（Dを右）</h3><pre>c1 —β1→ c2 —α3→ c3 —β2→ c4
   —α2→ c5 —β3→ c6 —α1→ c1</pre><p class="note">数学的な正の境界向き（Dを左）とは逆。境界だけを逆向きに読む際に、曲線の保存された向きやεを勝手に反転しない。</p></div><div><table><thead><tr><th>点</th><th>α</th><th>β</th><th>ε</th><th>面</th><th>役割</th></tr></thead><tbody>'''+xsrows+'''</tbody></table></div></div></section>
<section class="card"><h2>Koda変換：Dの6辺を除き、6個の角を平滑化</h2><p>α上の通過を負、β上の通過を正とする。頂点タイプl/rは通過符号とは別で、この規約ではσ<sub>o</sub>=−ε、l=−1、r=+1。</p><table><thead><tr><th>通る曲線</th><th>Dの外の経路</th><th>E-dataに追加する語</th></tr></thead><tbody>'''+blocksrows+'''</tbody></table><p>'''+wordhtml(actual)+'''</p><div class="two"><div><h3>o-graphのGauss表示</h3>'''+''.join(g)+'''<p class="note">円周の順序と頂点対の同定を示す図。弦同士の交差は新しいo-graph頂点ではない。実装用には10本の有向辺と元Heegaard辺への対応をKoda JSONに保存した。</p></div><div><h3>一致しない理由</h3><p>計算結果には正の連続ブロック <code>+3,+1,+4</code> がある。提示値には長さ3の同符号ブロックが一つもない。この性質は巡回移動・番号変更・全体の走査逆転・全通過符号反転でも解消しない。</p><table><thead><tr><th></th><th>正ブロック長</th><th>負ブロック長</th></tr></thead><tbody><tr><th>画像から計算</th><td>3,1,1</td><td>2,2,1</td></tr><tr><th>提示値</th><td>2,2,1</td><td>2,2,1</td></tr></tbody></table><p class="note">これは同じ符号付きo-dataかの照合。o-graph moveを経た同値や、表す3次元多様体の非同相性を主張するものではない。</p></div></div></section>
<section class="card"><h2>前後の経路を落とさない保存形式</h2><p>外周Oと3個の穴の輪郭L,M,Rを補助曲線として保存した。輪郭で切ると、前面・背面のそれぞれが種数0・境界4成分の曲面になる。14個の前後遷移点で貼り合わせるため、「穴のどちら側を回るか」をレイアウトに決めさせない。</p><p>全体は25頂点・50辺・21円板面のcellulation。これを46個の三角形分割された円板patchに展開して、前回のmarked validatorへ入力した。補助輪郭を描画対象から外したHeegaard版でも、同じ基準曲面と細分を保持する。</p><table><thead><tr><th>検査</th><th>結果</th></tr></thead><tbody><tr><th>曲面</th><td>閉・連結・向き付き、χ=−4、g=3</td></tr><tr><th>αで切断</th><td>連結、g=0、境界6</td></tr><tr><th>βで切断</th><td>連結、g=0、境界6</td></tr><tr><th>補領域</th><td>開円板7個</td></tr><tr><th>D</th><td>閉包も埋込み円板、6曲線との連結arc条件・方向条件に適合</td></tr><tr><th>描画</th><td>有理座標による全三角形・全辺検査、小数化後の再検査に合格</td></tr><tr><th>再現試験</th><td>Python 25項目合格。Julia用入力と試験を同梱するが、この環境では未実行。</td></tr></tbody></table><details><summary>機械可読な小型coreを表示（markingは別JSONに保存）</summary><pre>'''+html.escape(json.dumps(core,ensure_ascii=False,indent=2))+'''</pre></details></section>
<footer>文献：Y. Koda, “Branched spines and Heegaard genus of 3-manifolds”, manuscripta mathematica 123 (2007), 285–299, DOI:10.1007/s00229-007-0097-z。p.292: joining diskの定義、p.293: Fig.6のE-dataとα負・β正の規約。元画像はユーザー提供。手作業の画像転写と、それに対する厳密検証を区別した。自動画像認識や既存アプリへのインストールを行ったとはしていない。</footer></main>
<script>
const select=document.getElementById('choose');
function refresh(){const v=select.value;document.getElementById('original').style.opacity=v==='none'?'1':'.18';document.getElementById('labels').style.display=document.getElementById('showlabels').checked?'':'none';for(const p of document.querySelectorAll('.trace')){const show=(v==='all'||v===p.dataset.curve)&&(p.dataset.sheet==='F'?document.getElementById('front').checked:document.getElementById('back').checked);p.style.display=show?'':'none';}}
for(const x of document.querySelectorAll('.controls select,.controls input'))x.addEventListener('change',refresh);
for(const x of document.querySelectorAll('.point'))x.addEventListener('click',()=>{for(const r of document.querySelectorAll('.highlight'))r.classList.remove('highlight');const row=document.getElementById('row-'+x.dataset.point);row.classList.add('highlight');row.scrollIntoView({behavior:'smooth',block:'center'});});
refresh();
</script></html>'''
(root/'review.html').write_text(page,encoding='utf-8')
(root/'artifacts/gauss_diagram.svg').write_text(''.join(g),encoding='utf-8')
# also export a standalone actual o-graph with incidence, not just a word
og={'vertices':koda['vertices'],'directed_edges':koda['edges'],'gauss_components':koda['o_data'][0],
    'local_type_encoding':{'l':-1,'r':1},'meaning':'oriented typed normal o-graph; virtual projection crossings are not vertices'}
for v in og['vertices']:
    n=v['id'];v['branches']={}
    for visit,label in ((-n,'alpha'),(n,'beta')):
        j=actual.index(visit);v['branches'][label]={'visit':visit,'incoming':f'e{(j-1)%10+1}','outgoing':f'e{j+1}'}
B.dump(root/'data/figure6.ograph.json',og)
print('Wrote self-contained review.html and explicit o-graph incidence data.')
