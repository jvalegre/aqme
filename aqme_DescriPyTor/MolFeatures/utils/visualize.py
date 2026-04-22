import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
import sys
import pandas as pd
from typing import *
from enum import Enum
import os
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from utils  import help_functions as hf
except Exception as e :
    from .utils import help_functions as hf
# Now you can import from the parent directory



class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }

def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list

import re
from itertools import combinations

def plot_interactions(xyz_df, color, dipole_df=None, origin=None, sterimol_params=None):
    """
    Creates a 3D Plotly figure of the molecule (atoms + bonds) with optional dipole and Sterimol arrows.
    - Locks aspect to 'data' (no fake warping of small z).
    - Uses one coordinate array for all traces.
    - Correct visibility masks for toggles.
    """
    # ---------- helpers ----------
    def _parse_element(label: str) -> str:
        m = re.match(r"[A-Za-z]+", str(label).strip())
        if not m:
            s = str(label).strip()
        else:
            s = m.group(0)
        return s[0].upper() + s[1:].lower() if s else s

    def _resolve_origin(origin_arg, coords_arr):
        """
        origin_arg can be:
          - None: centroid of coords_arr
          - array-like shape (3,): xyz position
          - int: atom index (0-based)
          - list/array of ints: centroid of those atoms
        """
        if origin_arg is None:
            return coords_arr.mean(axis=0)
        ori = np.asarray(origin_arg)
        if ori.ndim == 1 and ori.shape[0] == 3 and np.isfinite(ori).all():
            return ori.astype(float)
        if np.issubdtype(ori.dtype, np.integer) and ori.ndim == 0:
            return coords_arr[int(ori)]
        if np.issubdtype(ori.dtype, np.integer) and ori.ndim == 1:
            return coords_arr[ori].mean(axis=0)
        raise ValueError("Invalid 'origin' argument. Use None, (3,), int, or list of ints.")

    def _planarity_metrics(coords_arr):
        """Optional diagnostic—unused in layout, handy to print if needed."""
        c = coords_arr.mean(0)
        X = coords_arr - c
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        n = Vt[-1]
        d = X @ n
        return dict(rms=float(np.sqrt((d**2).mean())), max_abs=float(np.max(np.abs(d))))

    # ---------- constants / lookups ----------
    try:
        atomic_radii_raw = GeneralConstants.COVALENT_RADII.value  # user env
        try:
            atomic_radii = dict(atomic_radii_raw)
        except Exception:
            atomic_radii = atomic_radii_raw
    except Exception:
        # minimal fallback
        atomic_radii = {"H": 0.31, "B": 0.85, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
                        "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
                        "Pd": 1.39, "Co": 1.26, "Ni": 1.24, "Fe": 1.24, "Cu": 1.32, "Zn": 1.22,
                        "Ag": 1.45, "Au": 1.36}

    cpk_colors = dict(
        C='black', F='green', H='white', N='blue', O='red', P='orange',
        S='yellow', Cl='green', Br='brown', I='purple',
        Ni='blue', Fe='red', Cu='orange', Zn='yellow', Ag='grey',
        Au='gold', Si='grey', B='pink', Pd='green', Co='pink'
    )

    # ---------- coordinates / atoms ----------
    coords = xyz_df[['x', 'y', 'z']].to_numpy(dtype=float)
    atoms_raw = xyz_df['atom'].astype(str).tolist()
    elements = [_parse_element(a) for a in atoms_raw]

    n = coords.shape[0]
    assert n == len(elements), f"Row count mismatch: coords={n}, atoms={len(elements)}"
    if not np.isfinite(coords).all():
        bad = np.argwhere(~np.isfinite(coords)).ravel().tolist()
        raise ValueError(f"Non-finite coordinates at rows: {bad}")

    # radii (safe defaults)
    default_r = 0.77
    radii = np.array([atomic_radii.get(el, default_r) for el in elements], dtype=float)

    # ---------- bonds (robust O(n^2)) ----------
    def get_bonds(thresh_scale=1.30, min_dist=0.10):
        bonds_local = {}
        for i, j in combinations(range(n), 2):
            dij = np.linalg.norm(coords[i] - coords[j])
            if dij <= min_dist:
                continue
            cutoff = (radii[i] + radii[j]) * thresh_scale
            if dij < cutoff:
                bonds_local[(i, j)] = round(dij, 2)
        return bonds_local

    bonds = get_bonds()

    # ---------- base traces ----------
    atom_colors = [cpk_colors.get(el, 'gray') for el in elements]
    hovertext = [f"{el} ({x:.3f},{y:.3f},{z:.3f})"
                 for el, (x, y, z) in zip(elements, coords)]
    atom_scatter = go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        mode='markers',
        marker=dict(color=atom_colors, size=5, line=dict(color='lightgray', width=2)),
        text=hovertext, name='atoms', hoverinfo='text'
    )

    bx, by, bz = [], [], []
    for (i, j) in bonds:
        xi, yi, zi = coords[i]
        xj, yj, zj = coords[j]
        bx += [xi, xj, None]
        by += [yi, yj, None]
        bz += [zi, zj, None]
    bond_trace = go.Scatter3d(
        x=bx, y=by, z=bz, mode='lines',
        line=dict(color=color, width=3), hoverinfo='none',
        name='bonds'
    )

    # ---------- annotation sets ----------
    annotations_idx = [
        dict(text=str(i + 1), x=coords[i, 0], y=coords[i, 1], z=coords[i, 2],
             showarrow=False, yshift=15, font=dict(color="blue"))
        for i in range(n)
    ]
    annotations_len = [
        dict(text=str(dist),
             x=(coords[i, 0] + coords[j, 0]) / 2,
             y=(coords[i, 1] + coords[j, 1]) / 2,
             z=(coords[i, 2] + coords[j, 2]) / 2,
             showarrow=False, yshift=10)
        for (i, j), dist in bonds.items()
    ]

    def add_traces():
        traces = []
        coords = xyz_df[['x','y','z']].to_numpy(float)
        center = coords.mean(axis=0)
        span = np.ptp(coords, axis=0)
        pad = np.median(span)*0.5   # how far beyond molecule you want to expand

        # create invisible boundary points
        x_dummy = [center[0] - pad, center[0] + pad]
        y_dummy = [center[1] - pad, center[1] + pad]
        z_dummy = [center[2] - pad, center[2] + pad]

        invisible_trace = go.Scatter3d(
            x=x_dummy,
            y=y_dummy,
            z=z_dummy,
            mode="markers",
            marker=dict(size=0.1, opacity=0),  # completely invisible
            hoverinfo="none",
            showlegend=False,
            name="frame_expander",
        )
        traces.append(invisible_trace)
            
        return traces
    
    # ---------- dipole arrows ----------
    def dipole_traces(dip_df, origin_arg, coords_for_scale, show_components=True):
        traces = []
        if dip_df is None or len(dip_df) == 0:
            return traces
        row = dip_df.iloc[0]
        vec = np.array([
            row.get("dipole_x", row.iloc[0]),
            row.get("dipole_y", row.iloc[1]),
            row.get("dipole_z", row.iloc[2]),
        ], dtype=float)
        if not np.all(np.isfinite(vec)):
            return traces

        tail = _resolve_origin(origin_arg, coords_for_scale)

        # auto-scale w.r.t. molecular span
        span = float(np.ptp(coords_for_scale, axis=0).max())
        dip_mag = float(np.linalg.norm(vec))
        scale = ((span / dip_mag) * 0.3 if dip_mag > 0 else 1.0) * 2.0

        comps = []
        if show_components:
            comps.extend([
                (np.array([vec[0] * scale, 0.0, 0.0]), "red",   "Dipole X"),
                (np.array([0.0, vec[1] * scale, 0.0]), "green", "Dipole Y"),
                (np.array([0.0, 0.0, vec[2] * scale]), "blue",  "Dipole Z"),
            ])
        comps.append((vec * scale, "purple", "Total Dipole"))

        for v, col, label in comps:
            L = float(np.linalg.norm(v))
            if L < 1e-10:
                continue
            end = tail + v
            shaft_end = end - 0.15 * (end - tail) / L

            traces.append(go.Scatter3d(
                x=[tail[0], shaft_end[0]],
                y=[tail[1], shaft_end[1]],
                z=[tail[2], shaft_end[2]],
                mode="lines",
                line=dict(color=col, width=4),
                name=label,
                showlegend=True,
            ))
            traces.append(go.Cone(
                x=[end[0]], y=[end[1]], z=[end[2]],
                u=[v[0]], v=[v[1]], w=[v[2]],
                anchor="tip",
                sizemode="scaled",
                sizeref=0.2,
                showscale=False,
                colorscale=[[0, col], [1, col]],
                name=label,
                showlegend=False,
            ))

        # mark origin
        traces.append(go.Scatter3d(
            x=[tail[0]], y=[tail[1]], z=[tail[2]],
            mode="markers",
            marker=dict(size=5, color="black", symbol="circle"),
            name="Origin (tail)"
        ))
        return traces

    # ---------- Sterimol arrows ----------
    def sterimol_traces(params, origin_arg):
        """
        Expects 'B1_coords','B5_coords','L_coords' as 2D vectors (x,y) in current global axes.
        Anchors at resolved 3D origin; lifts vectors with z=0 in global frame.
        """
        if params is None:
            return []
        required = ['B1_coords', 'B5_coords', 'L_coords', 'B1_value', 'B5_value', 'L_value']
        if not all(k in params for k in required):
            return []

        try:
            tail = _resolve_origin(origin_arg, coords)
        except Exception:
            tail = coords.mean(axis=0)

        vecs = [
            (np.asarray(params['B1_coords'], float), 'forestgreen', params['B1_value'], 'B1'),
            (np.asarray(params['B5_coords'], float), 'firebrick',   params['B5_value'], 'B5'),
            (np.asarray(params['L_coords'],  float), 'steelblue',   params['L_value'], 'L'),
        ]

        traces = []
        for vec2d, col, mag, name in vecs:
            vec3d = np.array([vec2d[0], vec2d[1], 0.0], float)
            L = float(np.linalg.norm(vec3d))
            if L < 1e-8:
                continue
            end = tail + vec3d
            shaft_end = end - 0.1 * (vec3d / L)

            traces.append(go.Scatter3d(
                x=[tail[0], shaft_end[0]],
                y=[tail[1], shaft_end[1]],
                z=[tail[2], shaft_end[2]],
                mode='lines',
                line=dict(color=col, width=4),
                name=f"{name}-shaft",
                showlegend=True
            ))
            traces.append(go.Cone(
                x=[end[0]], y=[end[1]], z=[end[2]],
                u=[vec3d[0]], v=[vec3d[1]], w=[vec3d[2]],
                anchor='tip',
                sizemode='absolute',
                sizeref=L * 0.07,   # head ~7% of length
                showscale=False,
                colorscale=[[0, col], [1, col]],
                name=name
            ))
        return traces

    # ---------- assemble data ----------
    data = [bond_trace, atom_scatter]

    dip_trs = dipole_traces(dipole_df, origin, coords_for_scale=coords) if dipole_df is not None else []
    for tr in dip_trs:
        tr.visible = False
    data.extend(dip_trs)

    add_trace=add_traces()
    data.extend(add_trace)
    for tr in add_trace:
        tr.visible = True
    ster_trs = sterimol_traces(sterimol_params, origin) if sterimol_params is not None else []
    for tr in ster_trs:
        tr.visible = False
    data.extend(ster_trs)

    # ---------- visibility masks ----------
    n_base = 2
    n_dip = len(dip_trs)
    n_ster = len(ster_trs)
    n_total = n_base + n_dip + n_ster

    vis_base_only = [True]*n_base + [False]*n_dip + [False]*n_ster
    vis_dip_only  = [True]*n_base + [True ]*n_dip + [False]*n_ster
    vis_ster_only = [True]*n_base + [False]*n_dip + [True ]*n_ster
    vis_both      = [True]*n_base + [True ]*n_dip + [True ]*n_ster
    # make add_trace vis
    # vis_add_only  = [True]*n_base + [False]*n_dip + [False]*n_ster + [True ]*n_add
    # ---------- buttons ----------
    buttons = [
        dict(label='Atom indices', method='relayout', args=[{'scene.annotations': annotations_idx}]),
        dict(label='Bond lengths', method='relayout', args=[{'scene.annotations': annotations_len}]),
        dict(label='Both',         method='relayout', args=[{'scene.annotations': annotations_idx + annotations_len}]),
        dict(label='Hide annotations', method='relayout', args=[{'scene.annotations': []}]),
    ]

    if n_dip > 0 and n_ster > 0:
        buttons.extend([
            dict(label='Show dipole',   method='update', args=[{'visible': vis_dip_only},  {}]),
            dict(label='Show Sterimol', method='update', args=[{'visible': vis_ster_only}, {}]),
            dict(label='Show both',     method='update', args=[{'visible': vis_both},      {}]),
            dict(label='Hide arrows',   method='update', args=[{'visible': vis_base_only}, {}]),
        ])
    elif n_dip > 0:
        buttons.extend([
            dict(label='Show dipole', method='update', args=[{'visible': vis_dip_only},  {}]),
            dict(label='Hide dipole', method='update', args=[{'visible': vis_base_only}, {}]),
        ])
    elif n_ster > 0:
        buttons.extend([
            dict(label='Show Sterimol', method='update', args=[{'visible': vis_ster_only}, {}]),
            dict(label='Hide Sterimol', method='update', args=[{'visible': vis_base_only}, {}]),
        ])

    updatemenus = [dict(buttons=buttons, direction='down', xanchor='left', yanchor='top')]

    return data, annotations_idx, updatemenus









def choose_conformers_input():
    string=input('Enter the conformers numbers: ')
    conformer_numbers=string.split(' ')
    conformer_numbers=[int(i) for i in conformer_numbers]
    return conformer_numbers


def unite_buttons(buttons_list, ref_index=0):
    buttons_keys=buttons_list[0].keys()
    united_button=dict.fromkeys(buttons_keys)
    for key in buttons_keys:
        if key=='args':
            all_annotations=[buttons[key][0]['scene.annotations'] for buttons in buttons_list]
            united_annotations=list(zip(*all_annotations))
            united_button[key]=[{'scene.annotations': united_annotations}]
        else:
            united_button[key]=buttons_list[ref_index][key]
    return united_button

def unite_updatemenus(updatemenus_list, ref_index=0):
    menus_keys=updatemenus_list[ref_index][0].keys()
    united_updatemenus_list=dict.fromkeys(menus_keys)
    for key in menus_keys:
        if key=='buttons':
            buttons_list=[updatemenus[0].get(key) for updatemenus in updatemenus_list]
            buttons_num=len(buttons_list[0])   
            segregated_buttons=[]
            for i in range(buttons_num):
                type_buttons=[buttons[i] for buttons in buttons_list]
                segregated_buttons.append(type_buttons)
            buttons=[unite_buttons(buttons) for buttons in segregated_buttons]
            united_updatemenus_list[key]=buttons
        else:
            united_updatemenus_list[key]=updatemenus_list[ref_index][0][key]
    return [united_updatemenus_list]


def compare_molecules(coordinates_df_list: List[pd.DataFrame],conformer_numbers:List[int]=None):
    if conformer_numbers is None:
        conformer_numbers=choose_conformers_input()
    # Create a subplot with 3D scatter plot
    colors_list=['red','purple','blue','green','yellow','orange','brown','black','pink','cyan','magenta']
    new_coodinates_df_list=[coordinates_df_list[i] for i in conformer_numbers]
    new_coodinates_df_list=renumbering_df_list(new_coodinates_df_list)
    # coordinates_df_list=renumber_xyz_by_mcs(coordinates_df_list)  ##needs fixing , renumbering not working.
    xyz_df=(coordinates_df_list[conformer_numbers[0]])
    data_main, annotations_id_main, updatemenus = plot_interactions(xyz_df,'grey')
    updatemenus_list=[updatemenus]
    # Iterate through the conformer numbers and create traces for each conformer
    for  conformer_number,color in zip((conformer_numbers[1:]),colors_list):
        xyz_df = coordinates_df_list[conformer_number]
        data, annotations_id, updatemenus_main = plot_interactions(xyz_df,color)
        data_main += data
        annotations_id_main += annotations_id
        updatemenus_list.append(updatemenus_main)
    # Set axis parameters
    updatemenus = unite_updatemenus(updatemenus_list)

    coords = xyz_df[['x','y','z']].to_numpy(float)
    span = np.ptp(coords, axis=0)            # total extent along each axis
    center = coords.mean(axis=0)
    pad = max(span) * 0.8                    # extra space around molecule

    x_range = [center[0] - span[0]/2 - pad, center[0] + span[0]/2 + pad]
    y_range = [center[1] - span[1]/2 - pad, center[1] + span[1]/2 + pad]
    z_range = [center[2] - span[2]/2 - pad, center[2] + span[2]/2 + pad]

    axis_params = dict(
        showgrid=True,
        showbackground=False,
        showticklabels=True,
        zeroline=True,
        titlefont=dict(color='black', size=12),
    )

    # --- Layout with fixed large frame ---
    layout = dict(
        scene=dict(
            aspectmode='data',
            xaxis=dict(**axis_params, title="X Axis", range=x_range),
            yaxis=dict(**axis_params, title="Y Axis", range=y_range),
            zaxis=dict(**axis_params, title="Z Axis", range=z_range),
            annotations=annotations_id_main,
        ),
        scene_camera=dict(
            eye=dict(x=1.8, y=1.8, z=1.8),
            center=dict(x=0, y=0, z=0),
            projection=dict(type='orthographic')
        ),
        margin=dict(r=0, l=0, b=0, t=30),
        showlegend=True,
        updatemenus=updatemenus,
        uirevision="keep-camera"
    )

    fig = go.Figure(data=data_main, layout=layout)
    fig.show()
    return fig

import dash
from dash import html, dcc, Output, Input, State
import plotly.graph_objects as go
import pandas as pd

def show_single_molecule(molecule_name,xyz_df=None,dipole_df=None, origin=None,sterimol_params=None,color='black'):
    if xyz_df is None:
        xyz_df=hf.get_df_from_file(hf.choose_filename()[0])
    # Create a subplot with 3D scatter plot
    data_main, annotations_id_main, updatemenus = plot_interactions(xyz_df,color,dipole_df=dipole_df, origin=origin,sterimol_params=sterimol_params)

    axis_params = dict(
        showgrid=False,
        showbackground=False,
        showticklabels=False,
        zeroline=False
    )

    scene = dict(
        xaxis=dict(**axis_params, title=dict(text='X', font=dict(color='white'))),
        yaxis=dict(**axis_params, title=dict(text='Y', font=dict(color='white'))),
        zaxis=dict(**axis_params, title=dict(text='Z', font=dict(color='white'))),
        annotations=annotations_id_main
    )

    layout = dict(
        title=dict(
            text=molecule_name,
            x=0.5,
            y=0.9,
            xanchor='center',
            yanchor='top'
        ),
        scene=scene,
        margin=dict(r=0, l=0, b=0, t=0),
        showlegend=False,
        updatemenus=updatemenus
        )
    
    fig = go.Figure(data=data_main, layout=layout)
    html=fig.show()
    run_app(fig)

    
    return html


def run_app(figure):
        # Create a Dash app
    app = dash.Dash(__name__)

    # App layout
    app.layout = html.Div([
        dcc.Graph(id='molecule-plot', figure=figure), # Replace "Water" with your molecule
        html.Div(id='clicked-data', children=[]),
        html.Button('Save Clicked Atom', id='save-button', n_clicks=0),
        html.Div(id='saved-atoms', children=[])
    ])

    # Callback to display clicked data
    @app.callback(
        Output('clicked-data', 'children'),
        Input('molecule-plot', 'clickData'),
        prevent_initial_call=True
    )
    def display_click_data(clickData):
        if clickData:
            return f"Clicked Point: {clickData['points'][0]['pointIndex']}"
        return "Click on an atom."

    # Callback to save clicked atom index
    @app.callback(
        Output('saved-atoms', 'children'),
        Input('save-button', 'n_clicks'),
        State('molecule-plot', 'clickData'),
        State('saved-atoms', 'children'),
        prevent_initial_call=True
    )
    def save_clicked_atom(n_clicks, clickData, saved_atoms):
        if clickData:
            saved_atoms.append(clickData['points'][0]['pointIndex'])
        return f"Saved Atom Indices: {saved_atoms}"
        
import matplotlib.pyplot as plt





# Atom color map (CPK-like)
atom_colors = {
    'C': 'black', 'H': 'gray', 'O': 'red', 'N': 'blue', 'S': 'yellow',
    'Cl': 'green', 'F': 'green', 'Br': 'brown', 'I': 'purple', 'P': 'orange'
}


from matplotlib.patches import ConnectionPatch

def plot_b1_visualization(rotated_plane, edited_coordinates_df,
                          n_points=100, title="Rotated Plane Visualization"):
    """
    Improved visualization of B1/B5 analysis:
      • Thin, crisp circles colored by atom identity
      • Non-overlapping numeric labels with leader lines
      • Highlighted B1/B5 arrows and angle arc
    """
    atom_colors = {
    'C': 'black', 'H': 'gray', 'O': 'red', 'N': 'blue', 'S': 'yellow',
    'Cl': 'green', 'F': 'green', 'Br': 'brown', 'I': 'purple', 'P': 'orange'
    }
    # ------------------ Compute extremes ------------------
    max_x, min_x = np.max(rotated_plane[:, 0]), np.min(rotated_plane[:, 0])
    max_y, min_y = np.max(rotated_plane[:, 1]), np.min(rotated_plane[:, 1])
    avs = np.abs([max_x, min_x, max_y, min_y])
    min_val, min_index = np.min(avs), np.argmin(avs)

    # B1 arrow (the minimal extreme)
    b1_coords = np.array([
        (max_x, 0), (min_x, 0), (0, max_y), (0, min_y)
    ][min_index])

    # B5 arrow (farthest point)
    norms_sq = np.sum(rotated_plane**2, axis=1)
    b5_idx = np.argmax(norms_sq)
    b5_point = rotated_plane[b5_idx]
    b5_value = np.linalg.norm(b5_point)

    # Angles
    angle_b1 = np.arctan2(b1_coords[1], b1_coords[0]) % (2*np.pi)
    angle_b5 = np.arctan2(b5_point[1], b5_point[0]) % (2*np.pi)
    angle_diff = abs(angle_b5 - angle_b1)
    if angle_diff > np.pi:
        angle_diff = 2*np.pi - angle_diff
    angle_diff_deg = np.degrees(angle_diff)

    # ------------------ Figure ------------------
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # ------------------ Circles ------------------
    n_total = rotated_plane.shape[0]
    n_circles = n_total // n_points
    centers = [
        (idx, row['atom'], row['y'], row['z'], row['radius'])
        for idx, row in edited_coordinates_df.iterrows()
    ]

    for i in range(n_circles):
        pts = rotated_plane[i*n_points:(i+1)*n_points, :]
        closed = np.vstack([pts, pts[0]])

        atom_idx, atom_type, y0, z0, radius = centers[i]
        color = atom_colors.get(atom_type, 'black')

        # Thin circle outline
        ax.plot(closed[:, 0], closed[:, 1],
                color=color, lw=0.8, alpha=0.9, solid_joinstyle='round')

        # Label slightly outside circle with connector line
        mean_y, mean_z = pts.mean(axis=0)
        vec = np.array([mean_y, mean_z])
        norm = np.linalg.norm(vec)
        if norm == 0: norm = 1e-5
        offset = 0.25 * vec / norm
        label_pos = vec + offset

        ax.text(label_pos[0], label_pos[1], str(atom_idx),
                ha='center', va='center', fontsize=8, fontweight='bold',
                color='black', bbox=dict(boxstyle='circle,pad=0.2',
                                         fc='white', ec='none', alpha=0.6))
        # connector line
        ax.add_patch(ConnectionPatch(xyA=(mean_y, mean_z),
                                     xyB=(label_pos[0], label_pos[1]),
                                     coordsA='data', coordsB='data',
                                     arrowstyle='-',
                                     lw=0.4, color=color, alpha=0.6))

    # ------------------ Axes lines ------------------
    ax.axvline(x=max_x, color='darkred', ls='dashed', lw=0.8, alpha=0.4)
    ax.axvline(x=min_x, color='darkred', ls='dashed', lw=0.8, alpha=0.4)
    ax.axhline(y=max_y, color='darkgreen', ls='dashed', lw=0.8, alpha=0.4)
    ax.axhline(y=min_y, color='darkgreen', ls='dashed', lw=0.8, alpha=0.4)

    # ------------------ B1 & B5 arrows ------------------
    arrow_colors = ['black']*4
    arrow_colors[min_index] = '#00A36C'   # highlight B1
    # X extremes
    ax.arrow(0, 0, max_x, 0, head_width=0.08, color=arrow_colors[0], length_includes_head=True)
    ax.arrow(0, 0, min_x, 0, head_width=0.08, color=arrow_colors[1], length_includes_head=True)
    # Y extremes
    ax.arrow(0, 0, 0, max_y, head_width=0.08, color=arrow_colors[2], length_includes_head=True)
    ax.arrow(0, 0, 0, min_y, head_width=0.08, color=arrow_colors[3], length_includes_head=True)

    # B5 arrow
    ax.arrow(0, 0, b5_point[0], b5_point[1],
             head_width=0.08, color="#CD3333", length_includes_head=True, lw=1.2, alpha=0.9)

    # Labels
    ax.text(b1_coords[0]*0.55, b1_coords[1]*0.55, f"B1\n{min_val:.2f}",
            fontsize=10, ha='center', va='bottom', fontweight='bold', color='#00A36C')
    ax.text(b5_point[0]*0.65, b5_point[1]*0.65, f"B5\n{b5_value:.2f}",
            fontsize=10, ha='center', va='bottom', fontweight='bold', color='#CD3333')

    # ------------------ Angle arc ------------------
    arc = np.linspace(min(angle_b1, angle_b5), max(angle_b1, angle_b5), 100)
    arc_x, arc_y = 0.6*np.cos(arc), 0.6*np.sin(arc)
    ax.plot(arc_x, arc_y, color='gray', lw=1.0, alpha=0.7)
    mid_angle = (min(angle_b1, angle_b5) + max(angle_b1, angle_b5)) / 2
    ax.text(0.75*np.cos(mid_angle), 0.75*np.sin(mid_angle),
            f"{angle_diff_deg:.1f}°", fontsize=9, color='gray',
            ha='center', va='center', fontweight='bold')

    ax.grid(alpha=0.15)
    plt.show()



def generate_circle(center_x, center_y, radius, n_points=20):
    """
    Generate circle coordinates given a center and radius.
    Returns a DataFrame with columns 'x' and 'y'.
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    x = center_x + radius * np.cos(theta)
    y = center_y + radius * np.sin(theta)
    return np.column_stack((x, y))


def plot_L_B5_plane(edited_coordinates_df, sterimol_df, n_points=100, title="L–B1 Plane (Y–Z)"):
    """
    Plot Y–Z plane with B1, L, and B5 arrows.
    Circles colored by atom type; atom index labels near each circle (no box).
    Summary labels below the legend; include B1–B5 angle at bottom.
    """
    # Build Y–Z circles & centers
    circles, centers = [], []
    for atom_idx, r in edited_coordinates_df.iterrows():
        pts = generate_circle(r['y'], r['z'], r['radius'], n_points=n_points)
        circles.append(pts)
        centers.append((atom_idx, r['atom'], r['y'], r['z'], r['radius']))
    plane_yz = np.vstack(circles)

    # Extract sterimol values
    L_val   = float(sterimol_df['L'].iloc[0])
    B5_val  = float(sterimol_df['B5'].iloc[0])
    loc_B5  = float(sterimol_df['loc_B5'].iloc[0])
    angle   = float(sterimol_df['B1_B5_angle'].iloc[0]) if 'B1_B5_angle' in sterimol_df.columns else None

    # Compute extremes
    ys, zs = plane_yz[:,0], plane_yz[:,1]
    max_y, min_y = ys.max(), ys.min()
    max_z, min_z = zs.max(), zs.min()
    abs_ext = np.abs([max_y, min_y, max_z, min_z])

    # Compute B1
    i_b1 = int(np.argmin(abs_ext))
    if i_b1 == 0:
        b1_y, b1_z = max_y, 0.0
    elif i_b1 == 1:
        b1_y, b1_z = min_y, 0.0
    elif i_b1 == 2:
        b1_y, b1_z = 0.0, max_z
    else:
        b1_y, b1_z = 0.0, min_z
    B1_val = abs_ext[i_b1]

    # Plot setup
    fig, ax = plt.subplots(figsize=(8,8))
    colors = {}
    # Plot circles and atom labels
    for i, pts in enumerate(circles):
        closed = np.vstack([pts, pts[0]])
        atom_idx, atom_type, y0, z0, radius = centers[i]
        color = atom_colors.get(atom_type, 'black')
        colors[atom_type]=color
        ax.plot(closed[:,0], closed[:,1], color=color, lw=1.5)
        offset = 0
        ax.text(y0 + offset, z0 + offset, str(atom_idx),
                ha='center', va='center', fontsize=9, color='black')

    # Dashed extremes
    ax.axvline(max_y, color='steelblue', ls='--')
    ax.axvline(min_y, color='steelblue', ls='--')
    ax.axhline(max_z, color='firebrick', ls='--')
    ax.axhline(min_z, color='firebrick', ls='--')

    # Plot arrows
    ax.arrow(0, 0, b1_y, b1_z, head_width=0.1, length_includes_head=True,
             color='gold', lw=2)
    y_L = max_y if abs(max_y) >= abs(min_y) else min_y
    ax.arrow(0, 0, y_L, 0, head_width=0.1, length_includes_head=True,
             color='forestgreen', lw=2)
    row = edited_coordinates_df.iloc[(edited_coordinates_df['y'] - loc_B5).abs().argmin()]
    y5, z5 = row['y'], row['z']
    ax.arrow(0, 0, y5, z5, head_width=0.1, length_includes_head=True,
             color='firebrick', lw=2)

    # Legend for atom types (upper right)
    legend_handles = [Line2D([0], [0], color=color, lw=3, label=atom)
                      for atom, color in colors.items()]
    ax.legend(handles=legend_handles, title="Atom Types",
              loc='upper left', bbox_to_anchor=(1.02, 1))

    # Summary labels below legend
    start_y = 0.5  # adjust as needed
    dy      = 0.1
    for i, (label, color) in enumerate([
            (f"B1: {B1_val:.2f}", 'gold'),
            (f"L:  {L_val:.2f}", 'forestgreen'),
            (f"B5: {B5_val:.2f}", 'firebrick'),
        ]):
        ax.text(
            1.02,
            start_y - i*dy,
            label,
            transform=ax.transAxes,
            ha='left', va='bottom',
            color=color,
            fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, lw=1)
        )
   

    # Angle annotation at bottom center
    if angle is not None:
        ax.text(0.5, -0.1, f"B1–B5 angle: {angle:.1f}°",
                transform=ax.transAxes, ha='center', va='top',
                fontsize=12, fontweight='bold')

    ax.set_title(title)
    ax.set_xlabel("Y")
    ax.set_ylabel("Z")
    ax.set_aspect('equal', 'box')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    
    pass
