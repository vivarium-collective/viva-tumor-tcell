"""Interactive Plotly visualizations for the tumor-tcell showcase.

Two headline forms:
  * ``timeseries_figure`` — multi-series line chart (population / phenotype
    fractions over time across conditions).
  * ``spatial_animation_figure`` — animated scatter of cells (colored by
    type+state, sized by radius) over the IFNg field heatmap, with a play
    button + time slider — the interactive analogue of the paper's videos.

Colors: Okabe–Ito (colorblind-safe). Cell-state categories are ordered so the
weakest CVD-adjacent pair (pink↔green) is never neighboring in the legend, and
every marker carries a thin outline (contrast relief for the lighter hues).
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

SURFACE = '#fcfcfb'
INK = '#1b1b1a'
MUTED = '#6b6b68'
GRID = '#e7e7e3'

# cell_type + cell_state -> (label, color)
STATE_STYLE = {
    ('tumor', 'PDL1n'):     ('Tumor PDL1n (proliferative)', '#0072B2'),
    ('tumor', 'PDL1p'):     ('Tumor PDL1p (arrested)',      '#E69F00'),
    ('t-cell', 'PD1n'):     ('T cell PD1- (active)',        '#009E73'),
    ('t-cell', 'PD1p'):     ('T cell PD1+ (exhausted)',     '#D55E00'),
    ('dendritic', 'inactive'): ('Dendritic (inactive)',     '#CC79A7'),
    ('dendritic', 'active'):   ('Dendritic (active)',        '#56B4E9'),
}
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


def spatial_animation_figure(frames, bounds, title, field='IFNg', max_frames=48):
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
            colorscale='Purples', zmin=0, zmax=fmax, opacity=0.55,
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


def write_html(fig, path, title):
    fig.write_html(str(path), include_plotlyjs='cdn', full_html=True,
                   config={'displayModeBar': True, 'responsive': True},
                   default_width='100%', default_height='680px')
