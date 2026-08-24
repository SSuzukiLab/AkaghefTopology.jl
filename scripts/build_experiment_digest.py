#!/usr/bin/env python3
"""Concatenate the experiment layer into one HTML document.

    python3 julia/scripts/build_experiment_digest.py [out.html]

The experiment layer is the `.qmd` slot in `Execution/qmd思想.md`: the artifacts
that say *what was investigated with the library*, as opposed to `src/`, which
says what can be computed. That slot is currently occupied by `julia/examples/`
in the wrong substrate -- `.jl` rather than `.qmd` -- so this digest is what
there is to convert.

Everything here is derived, never hand-written:

  LIVE_SCRIPT_SCOPE.md   the 60-file inventory and its T/D/R/E classification
  LIVESCRIPT_PORTS.md    per-script verified / partial status and open work
  examples/*.jl          the ported artifacts themselves, concatenated verbatim
  artifacts/layout_check.json   figure validity, when it has been generated

Regenerate rather than edit. A transcribed copy goes stale the moment a port
lands, which is the failure `qmd思想.md` names.
"""

import html
import json
import os
import pathlib
import re
import subprocess
import sys
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPO = ROOT.parent


def read(path):
    p = ROOT / path
    return p.read_text(encoding="utf-8") if p.exists() else ""


def git(*args, default=""):
    try:
        return subprocess.run(["git", "-C", str(REPO), *args],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return default


def scope_rows():
    """(class, mlx path, evidence) for all 60 inventoried Live Scripts."""
    return re.findall(r"^\| ([TDRE]) \| `([^`]+)` \| (.+?) \|$", read("LIVE_SCRIPT_SCOPE.md"), re.M)


def port_rows():
    """mlx path -> (status, evidence, open work) for the ported scripts."""
    out = {}
    for line in read("LIVESCRIPT_PORTS.md").splitlines():
        m = re.match(r"^\| `([^`]+)` \| (verified|partial) \| (.*?) \| (.*?) \|$", line)
        if m:
            out[m.group(1)] = (m.group(2), m.group(3), m.group(4))
    return out


def entry_points(source):
    return re.findall(r"^function (run_[A-Za-z0-9_]*)", source, re.M)


def layout_status():
    path = ROOT / "artifacts" / "layout_check.json"
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    per = {}
    for item in data["instances"]:
        key = os.path.basename(item["group"])
        entry = per.setdefault(key, {"pass": 0, "fail": 0, "failing": []})
        if item["status"] == "pass":
            entry["pass"] += 1
        elif item["status"] != "skipped":
            entry["fail"] += 1
            entry["failing"].append(f'{item["name"]} ({",".join(item.get("failing") or [])})')
    return {"summary": (data["checked"], data["failing"]), "per": per}


E = html.escape


def code_block(text):
    return f'<pre><code>{E(text)}</code></pre>'


def build():
    rows = scope_rows()
    ports = port_rows()
    examples = {p.stem: p for p in sorted((ROOT / "examples").glob("*.jl"))}
    layout = layout_status()

    done, todo = [], []
    for cls, mlx, evidence in rows:
        stem = os.path.basename(mlx).replace(".mlx", "")
        (done if stem in examples else todo).append((cls, mlx, evidence, stem))

    matched = {stem for _, _, _, stem in done}
    support = [(stem, path) for stem, path in examples.items() if stem not in matched]

    parts = []
    add = parts.append

    total_lines = sum(len(examples[s].read_text(encoding="utf-8").splitlines()) for _, _, _, s in done)
    total_lines += sum(len(p.read_text(encoding="utf-8").splitlines()) for _, p in support)

    # ---- masthead -------------------------------------------------------
    add(f"""<header class="mast">
<div>
  <h1>実験層ダイジェスト</h1>
  <p>`.qmd` に移すべき実験成果物の全体。完了しているものは本文ごと、できていないものは棚卸しとして、
     機械的に一つの文書へ結合したもの。</p>
</div>
<div class="stamp">
  生成 <b>{datetime.now(timezone.utc).astimezone().strftime('%Y-%m-%d %H:%M')}</b><br>
  commit <b>{E(git('rev-parse', '--short', 'HEAD', default='?'))}</b><br>
  <b>{len(done)}</b> 完了 / <b>{len(todo)}</b> 未着手
</div>
</header>

<div class="glance">
  <div class="g"><span class="n n--ok">{len(done)}</span><span class="k">移植済み実験</span></div>
  <div class="g"><span class="n n--bad">{len(todo)}</span><span class="k">未移植</span></div>
  <div class="g"><span class="n n--hold">{len(support)}</span><span class="k">補助ファイル</span></div>
  <div class="g"><span class="n">{total_lines:,}</span><span class="k">結合行数</span></div>
  <div class="g"><span class="n n--unk">0</span><span class="k">.qmd 化済み</span></div>
</div>

<section>
  <h2 id="premise">この文書の位置づけ</h2>
  <p class="lede">`qmd思想.md` の縦分業では、要約 → 論文 → 行間埋めノート → <b>実験群</b> → 実装 のうち、
  実験群だけが executable であるべき層で、その標準形式が `.qmd` である。</p>
  <p>現状その層は <code>julia/examples/*.jl</code> が占めている。移植は <code>.mlx</code> → <code>.qmd</code> ではなく
  <code>.mlx</code> → <code>.jl</code> で行われたため、<b>問い・条件・観測・解釈が失われ、実行部だけが残っている</b>。
  下の Part I はその実行部の全文であり、<code>.qmd</code> へ変換するときの素材そのものになる。</p>
  <p><b>0 件が <code>.qmd</code> 化済み</b>である。Quarto はこのマシンに導入されていない。</p>
</section>
""")

    # ---- part I ---------------------------------------------------------
    verified = [d for d in done if ports.get(d[1], ("", "", ""))[0] == "verified"]
    partial = [d for d in done if ports.get(d[1], ("", "", ""))[0] == "partial"]

    add(f"""<section>
  <h2 id="done">Part I — 完了している実験成果物 <span class="count">{len(done)}</span></h2>
  <p class="lede">分類 T（topology ワークフロー）と D（それが依存する代数）の全件。
  <b>verified {len(verified)} / partial {len(partial)}</b>。この status は
  <code>LIVESCRIPT_PORTS.md</code> 由来で、今回検算していない。</p>
  <p class="toolbar"><button id="expand" type="button">すべての本文を開く</button></p>
""")

    for group_name, group in (("verified", verified), ("partial", partial)):
        if not group:
            continue
        add(f'<h3 class="grouphead">{group_name} <span class="count">{len(group)}</span></h3>')
        for cls, mlx, evidence, stem in group:
            source = examples[stem].read_text(encoding="utf-8")
            status, port_ev, open_work = ports.get(mlx, ("—", "", ""))
            eps = entry_points(source)
            figs = ""
            if layout:
                key = os.path.basename(os.path.dirname(mlx)) or stem
                for gk, val in layout["per"].items():
                    if gk in mlx or gk == stem:
                        figs = f'図 {val["pass"]} 通過 / {val["fail"]} 違反'
                        if val["failing"]:
                            figs += " — " + E(", ".join(val["failing"]))
                        break
            add(f"""<article class="item">
  <div class="head">
    <div class="path">{E(mlx)}</div>
    <div class="badges"><span class="chip c--{'ok' if status=='verified' else 'hold'}">{status}</span>
      <span class="chip c--none">{cls}</span></div>
  </div>
  <dl class="meta">
    <dt>Julia</dt><dd><code>examples/{E(stem)}.jl</code> — {len(source.splitlines())} 行</dd>
    <dt>entry</dt><dd>{' '.join(f'<code>{E(e)}()</code>' for e in eps) or '—'}</dd>
    <dt>根拠</dt><dd>{E(port_ev) or '—'}</dd>
    <dt>未完了</dt><dd>{E(open_work) or '—'}</dd>
    {f'<dt>図</dt><dd>{figs}</dd>' if figs else ''}
  </dl>
  <details><summary>本文 {len(source.splitlines())} 行</summary>{code_block(source)}</details>
</article>""")

    if support:
        add(f'<h3 class="grouphead">補助ファイル <span class="count">{len(support)}</span></h3>'
            '<p class="note">対応する <code>.mlx</code> を持たない。集約や再生用のため、'
            '実験成果物ではなく実行基盤に属する。</p>')
        for stem, path in support:
            source = path.read_text(encoding="utf-8")
            eps = entry_points(source)
            add(f"""<article class="item">
  <div class="head"><div class="path">examples/{E(stem)}.jl</div>
    <div class="badges"><span class="chip c--none">support</span></div></div>
  <dl class="meta">
    <dt>entry</dt><dd>{' '.join(f'<code>{E(e)}()</code>' for e in eps) or '—'}</dd>
  </dl>
  <details><summary>本文 {len(source.splitlines())} 行</summary>{code_block(source)}</details>
</article>""")

    add("</section>")

    # ---- part II --------------------------------------------------------
    by_dir = {}
    for cls, mlx, evidence, stem in todo:
        by_dir.setdefault(os.path.dirname(mlx) or ".", []).append((cls, mlx, evidence))

    add(f"""<section>
  <h2 id="todo">Part II — できていないもの <span class="count">{len(todo)}</span></h2>
  <p class="lede">Julia 側に対応物が存在しない <code>.mlx</code>。全件が分類 R（topology API を呼ばない研究）
  または E（空・生成物）で、<b>T と D の欠けはゼロ</b>。</p>
  <p>プラットフォーム軸では移植不要である（Sage を呼ばないので Windows MATLAB でそのまま動く）。
  substrate 軸では対象に残る（<code>.mlx</code> は WikiLink を張れず、差分が取れず、agent が読めない）。
  <b>不要ではなく、急がない。</b></p>
""")
    for directory in sorted(by_dir):
        items = sorted(by_dir[directory], key=lambda r: r[1])
        add(f'<h3 class="grouphead">{E(directory)} <span class="count">{len(items)}</span></h3>'
            '<div class="scroll"><table><thead><tr><th style="width:44px">分類</th>'
            '<th class="mono">Live Script</th><th>内容</th></tr></thead><tbody>')
        for cls, mlx, evidence in items:
            add(f'<tr><td class="mono">{cls}</td>'
                f'<td class="mono">{E(os.path.basename(mlx))}</td>'
                f'<td>{E(re.sub("`", "", evidence))}</td></tr>')
        add("</tbody></table></div>")
    add("</section>")

    # ---- part III -------------------------------------------------------
    figsum = ""
    if layout:
        checked, failing = layout["summary"]
        figsum = (f'図の性質検査は {checked} インスタンス中 {failing} 件が違反。'
                  if failing else f'図の性質検査は {checked} インスタンス全件が通過。')
    else:
        figsum = ('<code>artifacts/layout_check.json</code> が未生成のため図の状態は不明。'
                  '<code>scripts/check_all_layouts.jl</code> を回すこと。')

    add(f"""<section>
  <h2 id="structural">Part III — 層そのものが未達な点</h2>
  <p class="lede">個々の移植とは別に、実験層の形式として満たしていないもの。</p>
  <ul class="plain">
    <li><span class="tag">substrate</span><div><b>0 件が <code>.qmd</code> 化されていない。</b>
      Quarto が未導入で、問い・条件・観測・解釈を持つ実験成果物は一つも存在しない。
      Part I にあるのは実行部だけである。</div></li>
    <li><span class="tag">二重責務</span><div><b><code>examples/</code> が実験記録とテスト fixture を兼ねている。</b>
      <code>test/runtests.jl</code> が {len(set(sum((entry_points(p.read_text(encoding='utf-8')) for p in examples.values()), [])))} 個の
      <code>run_*()</code> を参照する。<code>.qmd</code> へ移すには入力データを宣言的に切り出し、
      テストをそちらへ張り替える順序が要る。</div></li>
    <li><span class="tag">置き場所</span><div><b><code>.qmd</code> の物理的な位置が未決定。</b>
      <code>Execution/</code> 直下か <code>topology/experiments/</code> か。
      公開リポジトリと非公開 vault にまたがる WikiLink の扱いと連動する。</div></li>
    <li><span class="tag">図</span><div>{figsum}</div></li>
    <li><span class="tag">正本</span><div><b><code>qmd思想.md</code> が git 管理下にない。</b>
      設計の正本が版管理されていない。</div></li>
  </ul>
</section>""")

    return "\n".join(parts)


STYLE = """
<style>
:root{--ground:#f7f7f5;--panel:#fff;--ink:#12161d;--muted:#626a78;--faint:#8d95a2;
--rule:#d5d8dd;--rule-strong:#b4b9c2;--accent:#1b3bd0;--accent-soft:#e5e9fa;
--ok:#1a7346;--ok-soft:#dcece3;--hold:#8a6410;--hold-soft:#f2e9d5;
--bad:#ab2a1e;--bad-soft:#f6e0dd;
--display:"Iowan Old Style","Palatino Linotype",Palatino,Georgia,serif;
--body:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
--data:ui-monospace,SFMono-Regular,Menlo,"Cascadia Mono",Consolas,monospace}
@media(prefers-color-scheme:dark){:root{--ground:#0e1116;--panel:#151a21;--ink:#e6e9ef;
--muted:#939cab;--faint:#6f7887;--rule:#262c36;--rule-strong:#39414e;--accent:#6f8cff;
--accent-soft:#1a2340;--ok:#4cb583;--ok-soft:#16301f;--hold:#d3a441;--hold-soft:#2e2716;
--bad:#e4726a;--bad-soft:#331b1a}}
:root[data-theme=dark]{--ground:#0e1116;--panel:#151a21;--ink:#e6e9ef;--muted:#939cab;
--faint:#6f7887;--rule:#262c36;--rule-strong:#39414e;--accent:#6f8cff;--accent-soft:#1a2340;
--ok:#4cb583;--ok-soft:#16301f;--hold:#d3a441;--hold-soft:#2e2716;--bad:#e4726a;--bad-soft:#331b1a}
:root[data-theme=light]{--ground:#f7f7f5;--panel:#fff;--ink:#12161d;--muted:#626a78;
--faint:#8d95a2;--rule:#d5d8dd;--rule-strong:#b4b9c2;--accent:#1b3bd0;--accent-soft:#e5e9fa;
--ok:#1a7346;--ok-soft:#dcece3;--hold:#8a6410;--hold-soft:#f2e9d5;--bad:#ab2a1e;--bad-soft:#f6e0dd}
body{background:var(--ground);color:var(--ink);font-family:var(--body);font-size:16px;
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1040px;margin:0 auto;padding:0 24px 96px}
.mast{padding:56px 0 28px;border-bottom:2px solid var(--ink);display:flex;flex-wrap:wrap;
gap:28px;align-items:flex-end;justify-content:space-between}
.mast h1{font-family:var(--display);font-size:clamp(30px,4.4vw,44px);line-height:1.1;
font-weight:600;letter-spacing:-.01em;margin:0 0 10px;text-wrap:balance}
.mast p{margin:0;color:var(--muted);max-width:56ch}
.stamp{font-family:var(--data);font-size:12px;line-height:1.7;color:var(--muted);
text-align:right;white-space:nowrap}
.stamp b{color:var(--ink);font-weight:600}
.glance{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));
border-bottom:1px solid var(--rule)}
.g{padding:22px 20px 20px 0;border-right:1px solid var(--rule)}
.g:last-child{border-right:0}
.g .n{font-family:var(--data);font-size:30px;font-weight:600;letter-spacing:-.02em;
font-variant-numeric:tabular-nums;display:block;line-height:1.1}
.g .k{font-family:var(--data);font-size:10.5px;letter-spacing:.09em;text-transform:uppercase;
color:var(--muted);display:block;margin-top:7px}
.n--ok{color:var(--ok)}.n--bad{color:var(--bad)}.n--hold{color:var(--hold)}.n--unk{color:var(--faint)}
section{padding-top:52px}
h2{font-family:var(--display);font-size:26px;font-weight:600;margin:0 0 8px;
letter-spacing:-.005em;text-wrap:balance}
h3.grouphead{font-family:var(--data);font-size:11px;letter-spacing:.11em;text-transform:uppercase;
color:var(--muted);font-weight:600;margin:38px 0 14px;padding-bottom:8px;
border-bottom:1px solid var(--rule-strong)}
.count{font-family:var(--data);font-size:.62em;color:var(--faint);font-weight:400;
letter-spacing:.04em;margin-left:8px}
h2 .count{font-size:.5em}
.lede{color:var(--muted);margin:0 0 20px;max-width:70ch}
p{max-width:70ch}
.note{font-size:13.5px;color:var(--muted);max-width:70ch}
code{font-family:var(--data);font-size:.88em;background:var(--accent-soft);
color:var(--accent);padding:1px 4px}
.item{border:1px solid var(--rule);background:var(--panel);margin-bottom:14px}
.head{display:flex;flex-wrap:wrap;gap:12px;align-items:baseline;justify-content:space-between;
padding:14px 18px;border-bottom:1px solid var(--rule)}
.path{font-family:var(--data);font-size:13px;font-weight:600;word-break:break-all}
.badges{display:flex;gap:6px;flex:none}
.chip{font-family:var(--data);font-size:10.5px;letter-spacing:.05em;text-transform:uppercase;
padding:3px 7px;white-space:nowrap;font-weight:600}
.c--ok{background:var(--ok-soft);color:var(--ok)}
.c--hold{background:var(--hold-soft);color:var(--hold)}
.c--none{background:var(--accent-soft);color:var(--accent)}
dl.meta{display:grid;grid-template-columns:74px 1fr;gap:8px 16px;margin:0;padding:14px 18px;
font-size:13.5px;line-height:1.55}
dl.meta dt{font-family:var(--data);font-size:10.5px;letter-spacing:.08em;text-transform:uppercase;
color:var(--muted);padding-top:3px}
dl.meta dd{margin:0;color:var(--ink)}
details{border-top:1px solid var(--rule)}
summary{cursor:pointer;padding:11px 18px;font-family:var(--data);font-size:11.5px;
letter-spacing:.05em;color:var(--muted);user-select:none}
summary:hover{color:var(--accent)}
summary:focus-visible{outline:2px solid var(--accent);outline-offset:-2px}
pre{margin:0;padding:16px 18px;overflow-x:auto;background:var(--ground);
border-top:1px solid var(--rule)}
pre code{background:none;color:var(--ink);padding:0;font-size:12.5px;line-height:1.6}
.toolbar{margin:0 0 18px}
button{font-family:var(--data);font-size:11.5px;letter-spacing:.05em;padding:7px 14px;
background:var(--panel);color:var(--accent);border:1px solid var(--accent);cursor:pointer}
button:hover{background:var(--accent-soft)}
button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.scroll{overflow-x:auto;border-top:1px solid var(--rule-strong);border-bottom:1px solid var(--rule-strong);
margin-bottom:8px}
table{width:100%;border-collapse:collapse;font-size:14px;min-width:560px}
th{text-align:left;font-family:var(--data);font-size:10.5px;letter-spacing:.09em;
text-transform:uppercase;color:var(--muted);font-weight:600;padding:11px 14px 11px 0;
border-bottom:1px solid var(--rule)}
td{padding:11px 14px 11px 0;border-bottom:1px solid var(--rule);vertical-align:top;line-height:1.5}
tr:last-child td{border-bottom:0}
td.mono,th.mono{font-family:var(--data);font-size:12.5px}
ul.plain{list-style:none;padding:0;margin:0;display:flex;flex-direction:column}
ul.plain li{padding:14px 0;border-bottom:1px solid var(--rule);display:grid;
grid-template-columns:96px 1fr;gap:18px;align-items:start;max-width:80ch}
ul.plain li:first-child{border-top:1px solid var(--rule)}
.tag{font-family:var(--data);font-size:11px;font-weight:600;color:var(--accent);padding-top:2px}
footer{margin-top:60px;padding-top:22px;border-top:2px solid var(--ink);
font-size:13px;color:var(--muted);font-family:var(--data);line-height:1.8}
footer p{max-width:78ch}
@media(max-width:620px){.g{border-right:0;border-bottom:1px solid var(--rule);padding-right:0}
.stamp{text-align:left}dl.meta{grid-template-columns:1fr;gap:2px}
ul.plain li{grid-template-columns:1fr;gap:4px}}
@media(prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style>
"""

SCRIPT = """
<script>
(function(){
  var btn=document.getElementById('expand');
  if(!btn)return;
  btn.addEventListener('click',function(){
    var all=document.querySelectorAll('details');
    var open=btn.dataset.state!=='open';
    all.forEach(function(d){d.open=open;});
    btn.dataset.state=open?'open':'closed';
    btn.textContent=open?'すべての本文を閉じる':'すべての本文を開く';
  });
})();
</script>
"""

FOOTER = """
<footer>
<p>この文書は <code>julia/scripts/build_experiment_digest.py</code> の生成物であり、手で書いた箇所は無い。
出典は <code>LIVE_SCRIPT_SCOPE.md</code>（60件の棚卸しと分類）、<code>LIVESCRIPT_PORTS.md</code>（status）、
<code>examples/*.jl</code>（本文）、<code>artifacts/layout_check.json</code>（図）。編集せず再生成すること。</p>
<p>status の verified / partial は前セッションの文書由来で、検算していない。</p>
</footer>
"""


def main(argv):
    body = build()
    out = pathlib.Path(argv[0]) if argv else (ROOT / "artifacts" / "experiment_digest.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("<title>実験層ダイジェスト</title>\n" + STYLE +
                   '<div class="wrap">\n' + body + FOOTER + "</div>\n" + SCRIPT,
                   encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
