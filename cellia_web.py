from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
from typing import Optional, Dict, Any
import textwrap
import os
import json


cap_template = go.layout.Template(
    layout=go.Layout(
        font=dict(family="Inter, 'Helvetica Neue', Arial, sans-serif", size=14, color="#1f2937"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#FFFFFF",
        hoverlabel=dict(bgcolor="#FFFFFF", font_size=13, font_family="Inter, Arial", bordercolor="rgba(0,0,0,0.1)", font=dict(color="#1f2937")),
        legend=dict(
            orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
            bgcolor="rgba(255,255,255,0.8)", bordercolor="rgba(0,0,0,0.08)", borderwidth=1
        ),
        xaxis=dict(showgrid=True, gridcolor="#e5e7eb", zeroline=False),
        yaxis=dict(showgrid=True, gridcolor="#e5e7eb", zeroline=False),
        coloraxis_colorbar=dict(outlinewidth=0, ticks="outside", ticklen=4)
    )
)
pio.templates["cap_light_pro"] = cap_template
pio.templates.default = "cap_light_pro"

PASTEL = (
    px.colors.qualitative.Pastel + px.colors.qualitative.Set3 +
        px.colors.qualitative.D3 + px.colors.qualitative.T10
)  


"""PASTEL = None
def _ensure_plotly_theme():
    cap_template = go.layout.Template(
        layout=go.Layout(
            font=dict(family="Inter, 'Helvetica Neue', Arial, sans-serif", size=14, color="#1f2937"),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="#FFFFFF",
            hoverlabel=dict(bgcolor="#FFFFFF", font_size=13, font_family="Inter, Arial", bordercolor="rgba(0,0,0,0.1)", font=dict(color="#1f2937")),
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1,
                bgcolor="rgba(255,255,255,0.8)", bordercolor="rgba(0,0,0,0.08)", borderwidth=1
            ),
            xaxis=dict(showgrid=True, gridcolor="#e5e7eb", zeroline=False),
            yaxis=dict(showgrid=True, gridcolor="#e5e7eb", zeroline=False),
            coloraxis_colorbar=dict(outlinewidth=0, ticks="outside", ticklen=4)
        )
    )
    pio.templates["cap_light_pro"] = cap_template
    pio.templates.default = "cap_light_pro"

    global PASTEL
    if PASTEL is not None:
        return
    PASTEL = (
        px.colors.qualitative.Pastel + px.colors.qualitative.Set3 +
            px.colors.qualitative.D3 + px.colors.qualitative.T10
    )  """

def _discrete_map(series):
    cats = list(pd.Series(series).astype(str).unique())
    return {c: PASTEL[i % len(PASTEL)] for i, c in enumerate(sorted(cats))}


def _parse_marker_expl(value: str) -> dict:
    if not isinstance(value, str) or not value.strip(): return {}
    try:
        obj = json.loads(value)
        if isinstance(obj, dict): return {str(k): str(v) for k,v in obj.items()}
    except: pass
    expl = {}
    for p in value.replace(';',',').split(','):
        if ':' in p:
            g, r = p.split(':', 1)
            if g.strip(): expl[g.strip()] = r.strip()
    return expl


def build_marker_info_from_uns(adata, cluster_key="cluster", num_top_k: int = 15, rationale_json: Optional[Dict[str, Any]] = None):
    df_m = pd.DataFrame(adata.uns["marker_list"]).copy()
    df_g = pd.DataFrame(adata.uns.get("GPT_annotation_db", pd.DataFrame(columns=["cluster","LLM_annotation","markers","evidence_score","evidence_reason","marker_explanations"]))).copy()
    
    if "gene" in df_m.columns:
        df_m["gene"] = df_m["gene"].str.strip()

    df_m["cluster"] = df_m["cluster"].astype(str)
    if not df_g.empty:
        df_g["cluster"] = df_g["cluster"].astype(str)

    def _split_markers(s):
        if pd.isna(s): return []
        return [x.strip() for x in str(s).split(",") if x.strip()]

    if 'markers' in df_g.columns:
        df_g["markers_list"] = df_g["markers"].map(_split_markers)
    else:
        df_g["markers_list"] = [[] for _ in range(len(df_g))]


    db_by_cluster = (
        df_m.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
            .groupby("cluster")
            .apply(lambda sub: {
                "db_genes": sub["gene"].tolist(), "top15_db": sub["gene"].head(num_top_k).tolist(),
                "stats": sub[["gene","avg_log2FC","pct.1","pct.2","p_val_adj"]].to_dict("records")
            }).to_dict()
    )

    def process_llm_group(sub):
        first_row = sub.iloc[0]
        explanations = {}
        if 'marker_explanations' in first_row and pd.notna(first_row['marker_explanations']):
            explanations = _parse_marker_expl(first_row['marker_explanations'])
            
        return pd.Series({
            "cell_type": first_row.get("LLM_annotation", "Unknown"),
            "confidence": first_row.get("evidence_score"),
            "rationale": first_row.get("evidence_reason", ""),
            "llm_markers": first_row.get("markers_list", []),
            "marker_explanations": explanations
        })

    llm_by_cluster = {}
    if not df_g.empty:
        llm_by_cluster = (
            df_g.groupby("cluster")
                .apply(process_llm_group)
                .to_dict('index')
        )

    marker_info = {}
    clusters = sorted(set(df_m["cluster"]) | set(df_g.get("cluster", set())) | set(adata.obs[cluster_key].astype(str)))
    rationale_map = {str(k): v for k, v in rationale_json.items()} if rationale_json else {}

    for cid in clusters:
        llm = llm_by_cluster.get(cid, {})
        db = db_by_cluster.get(cid, {})
        rj = rationale_map.get(cid, {})
        
        cell_type = rj.get("cell_type", llm.get("cell_type", "Unknown"))
        confidence = rj.get("evidence_score", llm.get("confidence"))
        rationale = rj.get("evidence_reason", llm.get("rationale", ""))
        
        expl = {}
        llm_explained_genes = set()

        if rj.get("marker_explanations"):
            expl.update(rj["marker_explanations"])
            llm_explained_genes.update(rj["marker_explanations"].keys())
        elif llm.get("marker_explanations"):
            expl.update(llm["marker_explanations"])
            llm_explained_genes.update(llm["marker_explanations"].keys())
            for gene in db.get("top15_db", []):
                expl.setdefault(gene, "from DB (top 15 by avg_log2FC)")
        else:
            # Fallback if no 'marker_explanations' column exists.
            llm_genes = llm.get("llm_markers", [])
            llm_explained_genes.update(llm_genes)
            for gene in llm_genes: 
                expl[gene] = "LLM-selected marker"
            for gene in db.get("top15_db", []): 
                expl.setdefault(gene, "from DB (top 15 by avg_log2FC)")

        marker_info[cid] = {
            "cell_type": cell_type, "confidence": confidence, "marker_explanations": expl,
            "rationale": rationale, "db_stats": db.get("stats", []), 
            "llm_markers": llm.get("llm_markers", []),
            # [FIX] Store the specifically explained genes to use for highlighting.
            "llm_explained_genes": list(llm_explained_genes),
            "db_genes": db.get("db_genes", []), "top15_db": db.get("top15_db", []),
        }
    return marker_info

def make_umap_df(adata, marker_info, cluster_key="cluster"):
    df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'cell_id'}, inplace=True)

    df["cluster"] = adata.obs[cluster_key].astype(str).values
    df["cell_type"] = df["cluster"].map(lambda c: marker_info.get(c, {}).get("cell_type", "Unknown"))
    df["annotation"] = df["cluster"] + " (" + df["cell_type"] + ")"
    
    def get_markers_str(cluster):
        genes = list(marker_info.get(cluster, {}).get("marker_explanations", {}).keys())[:10]
        return "<br>".join(genes) if genes else "N/A"
    
    df["top_markers_str"] = df["cluster"].apply(get_markers_str)
    
    return df

def precompute_cluster_gene_stats(adata, marker_info, cluster_key="cluster"):
    genes_set = {g for v in marker_info.values() for g in v.get("llm_markers", []) + v.get("top15_db", [])}
    genes = [g for g in genes_set if g in adata.var_names]
    if not genes: return {"genes": [], "avg": [], "pct": []}
    X = adata[:, genes].X.toarray() if hasattr(adata[:, genes].X, "toarray") else adata[:, genes].X
    df_expr = pd.DataFrame(X, columns=genes, index=adata.obs_names)
    df_expr[cluster_key] = adata.obs[cluster_key].astype(str).values
    avg = df_expr.groupby(cluster_key)[genes].mean().reset_index().to_dict("records")
    pct = (df_expr[genes] > 0).groupby(df_expr[cluster_key]).mean().reset_index().to_dict("records")
    return {"genes": genes, "avg": avg, "pct": pct}

def _wrap_label(s: str, width: int = 20):
    return "<br>".join(textwrap.wrap(str(s), width=width))


def make_main_umap(umap_df, color_by="annotation", highlight_cluster_id=None):
    is_discrete = umap_df[color_by].dtype == 'object' or str(umap_df[color_by].dtype).startswith('category')
    
    fig = px.scatter(
        umap_df, x="UMAP1", y="UMAP2", color=color_by,
        render_mode="webgl",
        hover_name="cell_id",
        custom_data=["cluster", "cell_type", "top_markers_str"],
        color_discrete_map=_discrete_map(umap_df[color_by]) if is_discrete else None
    )

    fig.update_traces(
        marker=dict(size=7, opacity=0.85, line=dict(width=0)),
        hovertemplate=(
            "<b>Cell: %{hovertext}</b><br>"
            "Cluster: %{customdata[0]}<br>"
            "Annotation: %{customdata[1]}<br>"
            "<br><b>Top Markers:</b><br>"
            "%{customdata[2]}"
            "<extra></extra>"
        )
    )
    
    if highlight_cluster_id:
        sub = umap_df[umap_df["cluster"] == str(highlight_cluster_id)]
        if not sub.empty:
            fig.add_trace(go.Scattergl(
                x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
                marker=dict(
                    size=20,
                    color="#4F46E5",
                    opacity=0.4,
                    line=dict(width=0)
                ),
                name=f"Selected · {highlight_cluster_id}", hoverinfo="skip", showlegend=False
            ))
    fig.update_layout(
        title=None,
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="right", x=1),
        margin=dict(l=40, r=40, t=40, b=40)
    )
    return fig

def make_mini_umap(umap_df):
    mini = px.scatter(
        umap_df, x="UMAP1", y="UMAP2", color="cluster", custom_data=["cluster"],
        render_mode="webgl", color_discrete_map=_discrete_map(umap_df["cluster"])
    )
    mini.update_layout(
        margin=dict(l=0,r=0,t=0,b=0), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        showlegend=False
    )
    mini.update_traces(marker=dict(size=5, opacity=0.7), hovertemplate="Cluster %{customdata[0]}<extra></extra>")
    return mini

def make_explanation_card(cluster_id, marker_info, num_top_k: int = 15):
    info = marker_info.get(cluster_id, {})
    ct, conf, rationale = info.get("cell_type", "Unknown"), info.get("confidence"), info.get("rationale", "No rationale provided.")
    expl = info.get("marker_explanations", {})
    top15 = info.get("top15_db", [])
    
    # [FIX] Use the new, more specific list of explained genes for highlighting.
    llm_highlight = set(info.get("llm_explained_genes", []))

    chips = []
    for gene in top15:
        gene_str = str(gene)
        is_llm_gene = gene_str in llm_highlight
        icon = html.I(className="fas fa-brain", title="LLM-selected marker") if is_llm_gene else None
        chip = html.Span([gene_str, icon], className="chip")
        chips.append(chip)
        
    items = [html.Li([html.B(g), html.Span(f": {note}")], className="expl-item") for g, note in expl.items()]
    
    return html.Div([
        html.Div([html.Div("Selected Cluster Annotation", className="eyebrow"),
                  html.H3(f"Cluster {cluster_id} → {ct}", className="card-title"),
                  html.Div(f"Confidence Score: {conf:.2f}" if conf is not None else "Confidence: –", className="meta")]),
        html.Hr(),
        html.Div([html.Div("Rationale", className="eyebrow"), html.P(rationale, className="rationale-text")], className="section"),
        html.Hr(),
        html.Div([html.Div(f"Top {num_top_k} Markers", className="eyebrow"),
                  html.Div(chips or html.Span("No markers found.", className="muted"), className="chip-row")], className="section"),
        html.Hr(),
        html.Div([html.Div("Marker Explanations", className="eyebrow"),
                  html.Ul(items, className="expl-list") or html.Span("No specific explanations provided.", className="muted")], className="section"),
    ], className="card card--explain")

def make_dotplot_from_cache(cluster_id, marker_info, cache, topk=15):
    info = marker_info.get(cluster_id, {})
    llm, dbk = info.get("llm_markers", []), info.get("top15_db", [])
    
    all_genes_for_cluster = (llm + [g for g in dbk if g not in llm])
    genes_to_plot = [g for g in all_genes_for_cluster if g in cache.get('genes', [])][:topk]

    if not genes_to_plot:
        return go.Figure().update_layout(
            annotations=[dict(text="No marker data available.", showarrow=False, font=dict(size=14, color="#6c757d"))], 
            xaxis=dict(visible=False), yaxis=dict(visible=False), plot_bgcolor="#f8f9fa"
        )

    avg, pct = pd.DataFrame(cache["avg"]), pd.DataFrame(cache["pct"])
    avg_sub = avg.set_index("cluster")[genes_to_plot]
    pct_sub = pct.set_index("cluster")[genes_to_plot]
    labels = [f"{cid} ({marker_info.get(str(cid), {}).get('cell_type', 'Unknown')})" for cid in avg_sub.index]
    avg_sub.index = pct_sub.index = labels

    df_long = (avg_sub.reset_index().melt('index', var_name="Gene", value_name="AvgExpr")
               .merge(pct_sub.reset_index().melt('index', var_name="Gene", value_name="PctExpr"), on=["index", "Gene"]))
    
    orrd_scale = px.colors.sequential.OrRd
    custom_colorscale = [[0.0, "#ffffff"], [1e-9, orrd_scale[0]]]
    custom_colorscale.extend([[i / (len(orrd_scale)-1), color] for i, color in enumerate(orrd_scale)])


    fig = px.scatter(
        df_long, x="Gene", y="index", size="PctExpr", color="AvgExpr",
        color_continuous_scale=custom_colorscale,
        size_max=15,
        custom_data=["AvgExpr", "PctExpr"],
        title=None
    )
    
    fig.update_traces(
        hovertemplate=(
            "<b>Gene: %{x}</b><br>" +
            "Cluster Annotation: %{y}<br>" +
            "<br>" +
            "Avg. Expression: %{customdata[0]:.3f}<br>" +
            "Expressed in Cells: %{customdata[1]:.1%}" +
            "<extra></extra>"
        ),
        marker=dict(
            line=dict(width=1, color='#641a1a')
        )
    )

    fig.update_layout(
        xaxis_title=None, yaxis_title=None,
        margin=dict(l=40, r=40, t=40, b=120),
        xaxis=dict(tickangle=-45, automargin=True),
        yaxis=dict(automargin=True, tickfont=dict(size=10)),
        coloraxis_colorbar=dict(
            title="Avg.<br>Expr."
        )
    )
    return fig

def _load_rationale_file(path: str) -> Optional[dict]:
    if not path or not os.path.exists(path): return None
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext in ('.json', '.jsonl'):
            with open(path, 'r', encoding='utf-8') as f: data = json.load(f)
            if isinstance(data, list): df = pd.DataFrame(data)
            elif isinstance(data, dict): return {str(k): v for k, v in data.items()}
            else: return None
        else:
            df = pd.read_csv(path, sep='\t' if ext in ('.tsv', '.tab') else ',', encoding='utf-8')
    except Exception as e:
        print(f"Error loading rationale file {path}: {e}")
        return None

    if 'cluster' not in df.columns: return None
    df['cluster'] = df['cluster'].astype(str)
    out = {}
    for _, row in df.iterrows():
        cid = row['cluster']
        expl = _parse_marker_expl(row['marker_explanations']) if 'marker_explanations' in df.columns and pd.notna(row.get('marker_explanations')) else {}
        out[cid] = {
            'cell_type': row.get('cell_type'), 'evidence_score': row.get('evidence_score'),
            'evidence_reason': row.get('evidence_reason'), 'marker_explanations': expl
        }
    return out


def launch_cap_style_app(adata, port=8051, debug=True, num_top_k: int = 15, rationale_json_path: Optional[str] = None):
    #_ensure_plotly_theme()
    rj = _load_rationale_file(rationale_json_path)
    marker_info = build_marker_info_from_uns(adata, cluster_key="cluster", num_top_k=num_top_k, rationale_json=rj)
    umap_df = make_umap_df(adata, marker_info, cluster_key="cluster")
    cache_stats = precompute_cluster_gene_stats(adata, marker_info, cluster_key="cluster")

    app = Dash(__name__, external_stylesheets=['https://use.fontawesome.com/releases/v5.15.4/css/all.css'])
    app.title = "Cell Annotation Explorer"

    app.index_string = """
    <!DOCTYPE html>
    <html>
        <head>
            {%metas%}
            <title>Cell Annotation Explorer</title>
            {%favicon%}
            {%css%}
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
            <style>
                :root {
                    --bg-main: #f8f9fa; --bg-card: #ffffff; --bg-sidebar: #ffffff;
                    --ink-main: #212529; --ink-muted: #6c757d; --ink-subtle: #adb5bd; 
                    --accent-blue: #0d6efd; --border-color: #dee2e6;
                    --shadow-sm: 0 1px 2px 0 rgb(0 0 0 / 0.05);
                    --shadow-md: 0 4px 6px -1px rgb(0 0 0 / 0.08), 0 2px 4px -2px rgb(0 0 0 / 0.05);
                    --radius: 0.5rem; --font-main: 'Inter', sans-serif;
                }
                html, body { background: var(--bg-main); color: var(--ink-main); font-family: var(--font-main); height: 100vh; margin: 0; overflow: hidden; font-size: 15px; }
                * { box-sizing: border-box; }
                
                .app-container {
                    display: grid;
                    grid-template-columns: 340px 1fr;
                    height: 100vh;
                }
                .sidebar {
                    background: var(--bg-sidebar);
                    border-right: 1px solid var(--border-color);
                    padding: 24px;
                    display: flex;
                    flex-direction: column;
                    gap: 24px;
                    overflow-y: auto;
                }
                
                .sidebar-header h2 { font-size: 1.3rem; margin:0; }
                .sidebar-header p { color: var(--ink-muted); margin: 4px 0 0 0; font-size: 0.9rem;}
                
                .main-content { overflow-y: auto; padding: 32px; }
                .main-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }

                .card { 
                    background: var(--bg-card); border: 1px solid var(--border-color); border-radius: var(--radius); 
                    box-shadow: var(--shadow-sm);
                }
                
                .graph-card { padding: 20px; margin-bottom: 24px; display: flex; flex-direction: column;}
                .graph-title { font-size: 1rem; font-weight: 600; margin: 0 0 16px 8px; color: var(--ink-main); }
                .card--explain { padding: 24px 28px; }
                .card-title { margin: 0.2rem 0 0.5rem; font-size: 1.15rem; font-weight: 600; }
                .eyebrow { font-size: 0.7rem; color: var(--ink-muted); text-transform: uppercase; letter-spacing: .08em; font-weight: 600; }
                .meta { font-size: 0.9rem; color: var(--ink-muted); }
                hr { border: none; height: 1px; background-color: var(--border-color); margin: 20px 0; }
                
                .chip-row { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 12px; }
                .chip { display: inline-flex; align-items: center; gap: 8px; padding: 4px 10px; border-radius: 999px; background: #f1f3f5; font-size: 0.8rem; font-weight: 500; }
                .chip .fa-brain { color: var(--accent-blue); }
                
                .expl-list { margin: 12px 0 0; padding: 0 0 0 18px; list-style-type: '— '; }
                .expl-item { margin: 6px 0; line-height: 1.5; padding-left: 8px;}
                .rationale-text { line-height: 1.6; color: var(--ink-main); margin-top: 12px;}
                
                .control-group { display: flex; flex-direction: column; gap: 16px; }
                .control-group .Select-control, .control-group .rc-slider { width: 100%; }
                
                #mini-umap-container { height: 300px; padding: 8px; }

            </style>
        </head>
        <body>{%app_entry%}<footer>{%config%}{%scripts%}{%renderer%}</footer></body>
    </html>
    """

    # --- Sidebar Components ---
    sidebar = html.Div([
        html.Div([
            html.H2("Annotation Explorer"),
            html.P("Interactive single-cell data viewer.")
        ], className="sidebar-header"),
        
        html.Div([
            html.Div("Display Options", className="eyebrow"),
            html.Div(className="control-group", style={'marginTop': '12px'}, children=[
                dcc.Dropdown(
                    id="color-by", 
                    options=[{"label": x.replace("_", " ").title(), "value": x} for x in ["annotation", "cluster", "cell_type"]], 
                    value="annotation", 
                    clearable=False
                ),
                html.Div([
                    html.Label("Top K Genes in Dot Plot", style={'fontSize': '0.8rem', 'fontWeight': '500'}),
                    dcc.Slider(
                        id="topk", min=5, max=15, step=1, value=15, 
                        marks={5:"5", 10:"10", 15:"15"}, 
                        tooltip={"placement": "bottom"}
                    )
                ], style={'marginTop': '8px'})
            ])
        ], className="card", style={'padding': '16px'}),

        html.Div([
            html.Div("Navigation", className="eyebrow"),
            html.Div(dcc.Graph(id="mini-umap", figure=make_mini_umap(umap_df), config={'displayModeBar': False}), id="mini-umap-container", className="card", style={'marginTop': '12px'})
        ]),

    ], className="sidebar")

    # --- Main Content Components ---
    main_content = html.Div([
        html.Div([
            html.Div([
                html.H3("UMAP Visualization", className="graph-title"),
                html.Div(dcc.Graph(id="main-umap"), style={'height': '550px'})
            ], className="card graph-card"),
            html.Div([
                html.H3("Marker Gene Expression", className="graph-title"),
                html.Div(dcc.Graph(id="dot-plot"), style={'height': '550px'})
            ], className="card graph-card"),
        ], className="main-grid"),
        
        html.Div(id="marker-explanation-box", style={"marginTop": "24px"}),
    ], className="main-content")


    app.layout = html.Div([sidebar, main_content], className="app-container")

    # --- Callback Functions ---
    @app.callback(Output("main-umap", "figure"), [Input("color-by", "value"), Input("mini-umap", "clickData")])
    def update_main(color_by, clickData):
        highlight = str(clickData["points"][0]["customdata"][0]) if clickData else None
        return make_main_umap(umap_df, color_by=color_by, highlight_cluster_id=highlight)

    @app.callback([Output("marker-explanation-box", "children"), Output("dot-plot", "figure")], [Input("mini-umap", "clickData"), Input("topk", "value")])
    def update_details(clickData, topk):
        if not clickData:
            initial_card = html.Div([
                html.Div("Getting Started", className="eyebrow"),
                html.P("Click a cluster in the navigation UMAP on the left to view its detailed annotation, rationale, and marker genes.", style={"marginTop":"12px", "lineHeight": "1.6"})
            ], className="card card--explain")
            fig = go.Figure().update_layout(plot_bgcolor="#f8f9fa", annotations=[dict(text="Select a cluster to view dot plot", showarrow=False, font=dict(size=14, color="#6c757d"))])
            return initial_card, fig
        
        cluster_id = str(clickData["points"][0]["customdata"][0])
        card = make_explanation_card(cluster_id, marker_info, num_top_k=num_top_k)
        fig = make_dotplot_from_cache(cluster_id, marker_info, cache_stats, topk=topk)
        return card, fig

    app.run(debug=debug, port=port)