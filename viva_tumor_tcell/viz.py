"""Interactive Plotly visualizations for the tumor-tcell showcase.

Two headline forms:
  * ``timeseries_figure`` — multi-series line chart (population / phenotype
    fractions over time across conditions).
  * ``spatial_animation_figure`` — animated scatter of cells (colored by
    type+state, sized by radius) over the IFNg field heatmap, with a play
    button + time slider — the interactive analogue of the paper's videos.

Cell-state colors reproduce tumor-tcell's TAG_COLORS and the IFNg field uses its
'YlOrBr' colormap, so the spatial figures match the paper's published
snapshots/video. Condition/line charts use an Okabe–Ito set.
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

SURFACE = '#fcfcfb'
INK = '#1b1b1a'
MUTED = '#6b6b68'
GRID = '#e7e7e3'

# cell_type + cell_state -> (label, color).
# Colors are tumor-tcell's TAG_COLORS (matplotlib named colors), so the spatial
# figures match the paper's published snapshots/video exactly:
#   PDL1n=indianred, PDL1p=skyblue, PD1n=darkorange, PD1p=limegreen,
#   dendritic inactive=black, active=gray.
STATE_STYLE = {
    ('tumor', 'PDL1n'):     ('Tumor PDL1n (proliferative)', '#CD5C5C'),   # indianred
    ('tumor', 'PDL1p'):     ('Tumor PDL1p (arrested)',      '#87CEEB'),   # skyblue
    ('t-cell', 'PD1n'):     ('T cell PD1- (active)',        '#FF8C00'),   # darkorange
    ('t-cell', 'PD1p'):     ('T cell PD1+ (exhausted)',     '#32CD32'),   # limegreen
    ('dendritic', 'inactive'): ('Dendritic (inactive)',     '#000000'),   # black
    ('dendritic', 'active'):   ('Dendritic (active)',        '#808080'),   # gray
}
# IFNg field colormap — matches tumor-tcell's snapshots ('YlOrBr').
FIELD_COLORSCALE = 'YlOrBr'
CONDITION_COLORS = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7']


def _layout(fig, title, xaxis, yaxis):
    fig.update_layout(
        title=dict(text=title, font=dict(size=18, color=INK)),
        paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
        font=dict(family='system-ui, -apple-system, Segoe UI, sans-serif',
                  size=13, color=INK),
        margin=dict(l=64, r=24, t=56, b=52),
        legend=dict(bgcolor='rgba(0,0,0,0)', font=dict(size=12, color=INK)),
        hovermode='x unified',
    )
    fig.update_xaxes(title_text=xaxis, gridcolor=GRID, zeroline=False,
                     linecolor=GRID, ticks='outside', tickcolor=GRID, color=MUTED)
    fig.update_yaxes(title_text=yaxis, gridcolor=GRID, zeroline=False,
                     linecolor=GRID, ticks='outside', tickcolor=GRID, color=MUTED)
    return fig


def timeseries_figure(times_h, series, title, yaxis='count', xaxis='time (h)'):
    """series: {label: [y over time]} — a 2px line per series, unified hover."""
    fig = go.Figure()
    for i, (label, ys) in enumerate(series.items()):
        color = CONDITION_COLORS[i % len(CONDITION_COLORS)]
        fig.add_trace(go.Scatter(
            x=times_h, y=ys, name=label, mode='lines',
            line=dict(color=color, width=2.5),
            hovertemplate=f'<b>{label}</b>: %{{y:.3g}}<extra></extra>'))
    return _layout(fig, title, xaxis, yaxis)


def _hex_to_rgba(hexc, a):
    h = hexc.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f'rgba({r},{g},{b},{a})'


def timeseries_band_figure(times_h, series_stats, title, yaxis='count', xaxis='time (h)'):
    """Mean line + ±1 std band per condition across replicate seeds.

    series_stats: {label: (mean_array, std_array)}.
    """
    fig = go.Figure()
    for i, (label, (mean, std)) in enumerate(series_stats.items()):
        color = CONDITION_COLORS[i % len(CONDITION_COLORS)]
        mean = np.asarray(mean); std = np.asarray(std)
        upper = list(mean + std); lower = list(mean - std)
        # shaded ±std band
        fig.add_trace(go.Scatter(
            x=list(times_h) + list(times_h)[::-1], y=upper + lower[::-1],
            fill='toself', fillcolor=_hex_to_rgba(color, 0.15),
            line=dict(width=0), hoverinfo='skip', showlegend=False))
        fig.add_trace(go.Scatter(
            x=times_h, y=list(mean), name=label, mode='lines',
            line=dict(color=color, width=2.5),
            hovertemplate=f'<b>{label}</b>: %{{y:.3g}}<extra></extra>'))
    return _layout(fig, title, xaxis, yaxis)


def _time_hours(frames):
    return [f['time'] / 3600.0 for f in frames]


def spatial_animation_figure(frames, bounds, title, field='IFNg', max_frames=90):
    """Animated scatter of cells over the field heatmap.

    frames: list of snapshots (from run.snapshot_run). bounds: (bx, by) µm.
    """
    bx, by = float(bounds[0]), float(bounds[1])
    # subsample frames for a smooth, lightweight animation
    step = max(1, len(frames) // max_frames)
    idx = list(range(0, len(frames), step))
    if idx[-1] != len(frames) - 1:
        idx.append(len(frames) - 1)
    sframes = [frames[i] for i in idx]

    # shared field color range (robust upper bound)
    fmax = 1e-9
    for f in sframes:
        arr = np.asarray(f['fields'].get(field, [[0.0]]))
        if arr.size:
            fmax = max(fmax, float(np.percentile(arr, 99)))

    def _traces(frame):
        arr = np.asarray(frame['fields'].get(field, np.zeros((1, 1))), dtype=float)
        # field stored (nx, ny) -> heatmap wants z[row=y, col=x]
        z = arr.T
        nx, ny = arr.shape
        heat = go.Heatmap(
            z=z, x=np.linspace(0, bx, nx), y=np.linspace(0, by, ny),
            colorscale=FIELD_COLORSCALE, zmin=0, zmax=fmax, opacity=0.55,
            showscale=True, colorbar=dict(title=f'{field}<br>(ng/mL)', len=0.6, thickness=12),
            hoverinfo='skip')
        traces = [heat]
        # one scatter per state category (stable legend)
        for (ctype, cstate), (label, color) in STATE_STYLE.items():
            xs, ys, sizes, texts = [], [], [], []
            for cid, c in frame['cells'].items():
                if c.get('cell_type') == ctype and (c.get('cell_state') or '') == cstate:
                    loc = c.get('location', [0, 0])
                    xs.append(loc[0]); ys.append(loc[1])
                    sizes.append(max(6, float(c.get('radius', 5)) * 1.4))
                    texts.append(f"{cid}<br>{label}")
            traces.append(go.Scatter(
                x=xs, y=ys, mode='markers', name=label,
                marker=dict(color=color, size=sizes, line=dict(color='#2b2b28', width=1),
                            opacity=0.95),
                text=texts, hovertemplate='%{text}<extra></extra>',
                showlegend=True))
        return traces

    fig = go.Figure(data=_traces(sframes[0]))
    fig.frames = [go.Frame(data=_traces(f), name=f"{f['time']/3600.0:.1f}") for f in sframes]

    steps = [dict(method='animate', label=f"{f['time']/3600.0:.1f}",
                  args=[[f"{f['time']/3600.0:.1f}"],
                        dict(mode='immediate', frame=dict(duration=0, redraw=True),
                             transition=dict(duration=0))])
             for f in sframes]
    fig.update_layout(
        title=dict(text=title, font=dict(size=18, color=INK)),
        paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
        font=dict(family='system-ui, sans-serif', size=13, color=INK),
        margin=dict(l=52, r=24, t=56, b=64),
        legend=dict(bgcolor='rgba(252,252,251,0.7)', font=dict(size=11), y=0.99),
        xaxis=dict(title='x (µm)', range=[0, bx], gridcolor=GRID, color=MUTED,
                   constrain='domain'),
        yaxis=dict(title='y (µm)', range=[0, by], gridcolor=GRID, color=MUTED,
                   scaleanchor='x', scaleratio=1),
        updatemenus=[dict(type='buttons', showactive=False, x=0.02, y=1.08, xanchor='left',
                          buttons=[
                              dict(label='▶ Play', method='animate',
                                   args=[None, dict(frame=dict(duration=180, redraw=True),
                                                    fromcurrent=True,
                                                    transition=dict(duration=0))]),
                              dict(label='❚❚ Pause', method='animate',
                                   args=[[None], dict(mode='immediate',
                                                      frame=dict(duration=0, redraw=False))])])],
        sliders=[dict(active=0, x=0.12, len=0.85, xanchor='left',
                      currentvalue=dict(prefix='t = ', suffix=' h', font=dict(size=12)),
                      steps=steps)],
    )
    return fig


# ---- figures that mirror tumor-tcell's analysis plots ----

# per-state line colors (reuse the Okabe-Ito state palette)
_POP_COLORS = {
    'tumor_total': '#1b1b1a', 'tumor_PDL1n': '#0072B2', 'tumor_PDL1p': '#E69F00',
    'tcell_total': '#4b4b48', 'tcell_PD1n': '#009E73', 'tcell_PD1p': '#D55E00',
    'dendritic_total': '#6b6b68', 'dendritic_inactive': '#CC79A7', 'dendritic_active': '#56B4E9',
}
_POP_LABELS = {
    'tumor_total': 'Tumor total', 'tumor_PDL1n': 'Tumor PDL1n', 'tumor_PDL1p': 'Tumor PDL1p',
    'tcell_total': 'T cell total', 'tcell_PD1n': 'T cell PD1-', 'tcell_PD1p': 'T cell PD1+',
    'dendritic_total': 'Dendritic total', 'dendritic_inactive': 'Dendritic inactive',
    'dendritic_active': 'Dendritic active',
}


def population_group_figure(times_h, populations, keys, title, yaxis='cell count'):
    """Total + per-state counts over time (mirrors population_group_plot).

    populations: list of per-tick count dicts (from run.population_counts / analysis_run).
    keys: which count keys to draw, e.g. ['tumor_total','tumor_PDL1n','tumor_PDL1p'].
    """
    fig = go.Figure()
    for k in keys:
        ys = [p.get(k, 0) for p in populations]
        color = _POP_COLORS.get(k, INK)
        dash = 'solid' if k.endswith('total') else 'solid'
        width = 3 if k.endswith('total') else 2
        fig.add_trace(go.Scatter(
            x=times_h, y=ys, name=_POP_LABELS.get(k, k), mode='lines',
            line=dict(color=color, width=width, dash=dash),
            hovertemplate=f'<b>{_POP_LABELS.get(k, k)}</b>: %{{y}}<extra></extra>'))
    return _layout(fig, title, 'time (h)', yaxis)


def divisions_figure(times_h, div_series, title='Cumulative divisions', cell_types=('tumor', 't-cell')):
    """Cumulative divisions over time per cell type (mirrors division_plot)."""
    fig = go.Figure()
    colors = {'tumor': '#0072B2', 't-cell': '#009E73', 'dendritic': '#CC79A7'}
    for ct in cell_types:
        ys = [d.get(ct, 0) for d in div_series]
        fig.add_trace(go.Scatter(
            x=times_h, y=ys, name=f'{ct} divisions', mode='lines',
            line=dict(color=colors.get(ct, INK), width=2.5),
            hovertemplate=f'<b>{ct}</b>: %{{y}} divisions<extra></extra>'))
    return _layout(fig, title, 'time (h)', 'cumulative divisions')


# death-reason -> (label, color); mirrors the tumor-tcell death subtypes
_DEATH_STYLE = {
    'apoptosis':        ('Tumor apoptosis', '#E69F00'),
    'Tcell_death':      ('Tumor killed by T cell', '#D55E00'),
    'PD1n_apoptosis':   ('T cell PD1- apoptosis', '#009E73'),
    'PD1p_apoptosis':   ('T cell PD1+ apoptosis', '#56B4E9'),
    'PD1p_PDL1_death':  ('T cell PD1+ death at PDL1+', '#CC79A7'),
}


def deaths_figure(times_h, death_series, title='Cumulative deaths by type'):
    """Cumulative deaths by death-reason type over time (mirrors death_group_plot)."""
    reasons = sorted({r for d in death_series for r in d})
    fig = go.Figure()
    for i, r in enumerate(reasons):
        label, color = _DEATH_STYLE.get(r, (r, CONDITION_COLORS[i % len(CONDITION_COLORS)]))
        ys = [d.get(r, 0) for d in death_series]
        fig.add_trace(go.Scatter(
            x=times_h, y=ys, name=label, mode='lines',
            line=dict(color=color, width=2.5),
            hovertemplate=f'<b>{label}</b>: %{{y}}<extra></extra>'))
    if not reasons:
        fig.add_annotation(text='no deaths recorded in this run', showarrow=False,
                           xref='paper', yref='paper', x=0.5, y=0.5, font=dict(color=MUTED))
    return _layout(fig, title, 'time (h)', 'cumulative deaths')


def snapshots_panel_figure(snapshots, bounds, title, field='IFNg'):
    """Multi-timepoint spatial panels: cells (colored by type+state, sized by
    radius) over the field heatmap — a static analogue of tumor-tcell's
    plot_snapshots (one column per timepoint)."""
    from plotly.subplots import make_subplots
    bx, by = float(bounds[0]), float(bounds[1])
    n = len(snapshots)
    titles = [f"t = {s['time']/3600.0:.1f} h" for s in snapshots]
    fig = make_subplots(rows=1, cols=n, subplot_titles=titles,
                        horizontal_spacing=0.02, shared_yaxes=True)
    fmax = max([float(np.percentile(np.asarray(s['fields'].get(field, [[0.0]])), 99))
                for s in snapshots] + [1e-9])
    seen_legend = set()
    for col, s in enumerate(snapshots, start=1):
        arr = np.asarray(s['fields'].get(field, np.zeros((1, 1))), dtype=float)
        fig.add_trace(go.Heatmap(z=arr.T, x=np.linspace(0, bx, arr.shape[0]),
                                 y=np.linspace(0, by, arr.shape[1]),
                                 colorscale=FIELD_COLORSCALE, zmin=0, zmax=fmax, opacity=0.5,
                                 showscale=False, hoverinfo='skip'), row=1, col=col)
        for (ct, cs), (label, color) in STATE_STYLE.items():
            xs, ys, sz = [], [], []
            for c in s['cells'].values():
                if c.get('cell_type') == ct and (c.get('cell_state') or '') == cs:
                    loc = c.get('location', [0, 0]); xs.append(loc[0]); ys.append(loc[1])
                    sz.append(max(5, float(c.get('radius', 5)) * 1.1))
            show = label not in seen_legend
            if xs:
                seen_legend.add(label)
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode='markers', name=label, legendgroup=label, showlegend=show,
                marker=dict(color=color, size=sz, line=dict(color='#2b2b28', width=0.7)),
                hovertemplate=f'{label}<extra></extra>'), row=1, col=col)
        fig.update_xaxes(range=[0, bx], showgrid=False, row=1, col=col,
                         scaleanchor=f'y{col if col > 1 else ""}', constrain='domain')
        fig.update_yaxes(range=[0, by], showgrid=False, row=1, col=col)
    fig.update_layout(title=dict(text=title, font=dict(size=17, color=INK)),
                      paper_bgcolor=SURFACE, plot_bgcolor=SURFACE, height=380,
                      font=dict(size=12, color=INK),
                      legend=dict(orientation='h', y=-0.12, font=dict(size=11)))
    return fig


def cytotoxicity_figure(times_h, cyto_mean, cyto_sem, title='Cytotoxicity vs. control (mean ± SEM)'):
    """% cytotoxicity over time with a mean ± SEM band (mirrors cytotoxicity_rep_plot)."""
    fig = go.Figure()
    color = '#D55E00'
    upper = list(np.asarray(cyto_mean) + np.asarray(cyto_sem))
    lower = list(np.asarray(cyto_mean) - np.asarray(cyto_sem))
    fig.add_trace(go.Scatter(x=list(times_h) + list(times_h)[::-1], y=upper + lower[::-1],
                             fill='toself', fillcolor=_hex_to_rgba(color, 0.15),
                             line=dict(width=0), hoverinfo='skip', showlegend=False))
    fig.add_trace(go.Scatter(x=times_h, y=list(cyto_mean), name='cytotoxicity', mode='lines',
                             line=dict(color=color, width=2.5),
                             hovertemplate='%{y:.1f}%<extra></extra>'))
    return _layout(fig, title, 'time (h)', 'cytotoxicity (%)')


def spatial_gif_html(snapshots, bounds, title, field='IFNg', fps=10, max_frames=120):
    """Render an animated GIF of the spatial dynamics (cells over the IFNg field)
    with the paper's TAG_COLORS / YlOrBr, and return it wrapped in an HTML page
    (autoplaying, looping <img>) — the tumor-tcell 'video' analogue.

    Cells are drawn as true ``patches.Circle`` in µm data-coordinates on an
    equal-aspect axis, so a 15 µm tumor is visibly larger than a 7.5 µm T cell
    (unlike a pixel-sized scatter marker), and every retained frame is rendered
    so motion is smooth rather than teleporting. Uses matplotlib (Agg) like the
    original's snapshots/video.
    """
    import base64
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib import animation, cm
    from matplotlib.colors import Normalize

    bx, by = float(bounds[0]), float(bounds[1])
    step = max(1, len(snapshots) // max_frames)
    frames = snapshots[::step]
    fmax = max([float(np.percentile(np.asarray(f['fields'].get(field, [[0.0]])), 99))
                for f in frames] + [1e-9])

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    fig.patch.set_facecolor('white')
    # colorbar drawn once on its own axes (survives ax.clear() each frame)
    norm = Normalize(vmin=0.0, vmax=fmax)
    sm = cm.ScalarMappable(norm=norm, cmap='YlOrBr')
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(f'{field} (ng/mL)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    def draw(frame):
        ax.clear()
        arr = np.asarray(frame['fields'].get(field, np.zeros((1, 1))), dtype=float)
        ax.imshow(arr.T, origin='lower', extent=[0, bx, 0, by], cmap='YlOrBr',
                  norm=norm, alpha=0.85, aspect='equal', zorder=0)
        for c in frame['cells'].values():
            style = STATE_STYLE.get((c.get('cell_type'), c.get('cell_state') or ''))
            if not style:
                continue
            loc = c.get('location', (0, 0))
            r = float(c.get('radius', 5.0))
            ax.add_patch(patches.Circle((loc[0], loc[1]), r, facecolor=style[1],
                                        edgecolor='white', linewidth=0.5, zorder=2))
        ax.set_xlim(0, bx); ax.set_ylim(0, by)
        ax.set_aspect('equal')
        ax.set_title(f"{title}\nt = {frame['time']/3600.0:.1f} h", fontsize=10)
        ax.set_xlabel('x (µm)', fontsize=9); ax.set_ylabel('y (µm)', fontsize=9)
        ax.tick_params(labelsize=8)
        return []

    anim = animation.FuncAnimation(fig, draw, frames=frames, blit=False)
    import os
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix='.gif', delete=False)
    tmp.close()
    try:
        anim.save(tmp.name, writer=animation.PillowWriter(fps=fps))
        plt.close(fig)
        with open(tmp.name, 'rb') as fh:
            b64 = base64.b64encode(fh.read()).decode('ascii')
    finally:
        os.unlink(tmp.name)
    legend = ''.join(
        f'<span style="display:inline-block;margin:0 8px 0 0;font:12px system-ui">'
        f'<span style="display:inline-block;width:10px;height:10px;border-radius:50%;'
        f'background:{color};border:1px solid #999;vertical-align:middle"></span> {label}</span>'
        for (label, color) in STATE_STYLE.values())
    return (
        f'<!doctype html><html><body style="margin:0;background:{SURFACE};'
        f'text-align:center;font-family:system-ui">'
        f'<div style="padding:6px">{legend}</div>'
        f'<img alt="{title}" style="max-width:100%;height:auto" '
        f'src="data:image/gif;base64,{b64}"/></body></html>')


def write_html_str(html, path):
    """Write a pre-rendered HTML string (e.g. from spatial_gif_html)."""
    from pathlib import Path
    Path(path).write_text(html, encoding='utf-8')


def write_html(fig, path, title):
    fig.write_html(str(path), include_plotlyjs='cdn', full_html=True,
                   config={'displayModeBar': True, 'responsive': True},
                   default_width='100%', default_height='680px')
