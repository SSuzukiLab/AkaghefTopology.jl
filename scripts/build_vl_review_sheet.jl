using Topology, JSON

const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"examples","VLExample.jl"))

html_escape(value)=replace(string(value),
    "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\""=>"&quot;")

function review_sheet(;output_dir=raw"C:/Users/Akaghef/AppData/Local/Temp/vl-review")
    plots_dir=joinpath(output_dir,"plots")
    mkpath(plots_dir)
    calls=run_vl_plot_calls()
    all_examples=copy(calls)
    virtual=vl_virtual_gauss_examples()
    push!(all_examples,(name="poincare",link=virtual["poincare"],bending=nothing))

    records=Any[]
    for (index,call) in enumerate(all_examples)
        slug=replace(call.name,r"[^A-Za-z0-9_-]"=>"_")
        plot_path=joinpath(plots_dir,"$(lpad(index,2,'0'))_$(slug).svg")
        plot_status="not-run"
        plot_error=""
        diagnostics=Any[]
        try
            plot_svg(call.link,plot_path;bending_numbers=call.bending)
            plot_status="generated"
            diagnostics=check_layout(call.link;bending_numbers=call.bending)
        catch error
            plot_status="error"
            plot_error=sprint(showerror,error)
        end
        push!(records,Dict(
            "index"=>index,
            "name"=>call.name,
            "plot_status"=>plot_status,
            "plot_error"=>plot_error,
            "diagnostic_count"=>length(diagnostics),
            "error_count"=>count(diagnostic -> diagnostic.severity===:error,diagnostics),
            "failing_rules"=>string.(unique(diagnostic.rule for diagnostic in diagnostics if diagnostic.severity===:error)),
            "plot_file"=>"plots/$(lpad(index,2,'0'))_$(slug).svg",
        ))
    end

    generated=count(record -> record["plot_status"]=="generated",records)
    errors=sum(record["error_count"] for record in records)
    cards=join([
        begin
            status=record["plot_status"]=="generated" ? "generated" : "ERROR"
            diagnostics=record["error_count"]==0 ? "layout: OK" :
                "layout errors: $(record["error_count"]) ($(join(record["failing_rules"],", ")))"
            error_block=isempty(record["plot_error"]) ? "" :
                "<pre class=\"error\">$(html_escape(record["plot_error"]))</pre>"
            """
            <article class="card" data-name="$(html_escape(record["name"]))">
              <div class="figure">
                <h2>$(lpad(record["index"],2,'0')) · $(html_escape(record["name"]))</h2>
                <img src="$(record["plot_file"])" alt="$(html_escape(record["name"])) plot">
              </div>
              <div class="diagnostics">
                <strong class="$status">$status</strong>
                <div>$diagnostics</div>
                $error_block
              </div>
              <div class="review">
                <label>Review</label>
                <textarea data-review="$(html_escape(record["name"]))" placeholder="確認メモ、修正点、保留理由…"></textarea>
              </div>
            </article>
            """
        end for record in records],"\n")
    html="""
    <!doctype html>
    <html lang="ja">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width,initial-scale=1">
      <title>VirtualLink review sheet</title>
      <style>
        :root { color-scheme: light; font-family: system-ui, sans-serif; }
        body { margin: 0; padding: 24px; color: #172033; background: #f5f7fb; }
        header { position: sticky; top: 0; z-index: 2; padding: 16px 20px; margin: -24px -24px 20px;
                 background: rgba(245,247,251,.96); border-bottom: 1px solid #d7ddea; }
        h1 { margin: 0 0 6px; font-size: 24px; }
        .summary { color: #536078; }
        .grid { display: grid; gap: 16px; }
        .card { display: grid; grid-template-columns: minmax(360px, 1.5fr) 220px minmax(240px, .8fr);
                gap: 16px; align-items: stretch; padding: 14px; background: white;
                border: 1px solid #d7ddea; border-radius: 12px; box-shadow: 0 2px 8px #17203312; }
        h2 { margin: 0 0 8px; font-size: 16px; }
        .figure { min-width: 0; }
        .figure img { display: block; width: 100%; height: 300px; object-fit: contain; background: #fff; border: 1px solid #edf0f6; }
        .diagnostics { padding: 12px; background: #f7f9fd; border-radius: 8px; line-height: 1.6; }
        .generated { color: #087443; }
        .ERROR { color: #b42318; }
        .error { white-space: pre-wrap; overflow-wrap: anywhere; font-size: 11px; color: #b42318; }
        .review { display: flex; flex-direction: column; gap: 8px; }
        .review label { font-weight: 700; }
        textarea { flex: 1; min-height: 180px; resize: vertical; padding: 10px; border: 1px solid #c8d0df; border-radius: 8px; font: inherit; }
        @media (max-width: 900px) { .card { grid-template-columns: 1fr; } .figure img { height: 260px; } }
      </style>
    </head>
    <body>
      <header>
        <h1>VirtualLink review sheet</h1>
        <div class="summary">examples: $(length(records)) · SVG generated: $generated · layout diagnostic errors: $errors</div>
        <div class="summary">各カード：図 ｜ 幾何診断 ｜ Review欄</div>
      </header>
      <main class="grid">$cards</main>
      <script>
        for (const area of document.querySelectorAll('textarea[data-review]')) {
          const key = 'VirtualLinkReview:' + area.dataset.review;
          area.value = localStorage.getItem(key) || '';
          area.addEventListener('input', () => localStorage.setItem(key, area.value));
        }
      </script>
    </body>
    </html>
    """
    html_path=joinpath(output_dir,"review.html")
    write(html_path,html)
    json_path=joinpath(output_dir,"review-status.json")
    open(json_path,"w") do io
        JSON.print(io,Dict("examples"=>records,"generated"=>generated,"layout_errors"=>errors),2)
    end
    println("review_html=",html_path)
    println("status_json=",json_path)
    println("examples=",length(records)," generated=",generated," layout_errors=",errors)
    html_path
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    review_sheet()
end
