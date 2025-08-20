"""
Python implementation of R script for optimizing color palettes
Maximizes harmonic mean pairwise color distances in perceptually uniform colorspaces
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from scipy import optimize
import colorsys
from colorspacious import cspace_convert
from matplotlib.colors import to_rgb, to_hex, hex2color
import matplotlib.cm as cm

# Helper functions for coordinate transformations
def polar2cart(t, r):
    """Convert polar to Cartesian coordinates"""
    return np.array([r * np.cos(t), r * np.sin(t)])

def logit(p):
    """Logit transformation"""
    return np.log(p / (1 - p))

def invlogit(x):
    """Inverse logit transformation"""
    return np.exp(x) / (1 + np.exp(x))

def cart2polar(x, y):
    """Convert Cartesian to polar coordinates"""
    return np.array([np.arctan2(y, x), np.sqrt(x**2 + y**2)])

def shrink(x, a=0.1):
    """Shrink a range by a factor"""
    return np.array([x[0] + (x[1] - x[0])/2 * a, x[1] - (x[1] - x[0])/2 * a])

# Convert HSL to hex (matplotlib compatible)
def hsl_to_hex(h, s, l):
    """
    Convert HSL color to hex format
    h: hue (0-360)
    s: saturation (0-100%)
    l: lightness (0-100%)
    """
    # Convert to 0-1 range
    h = h / 360.0
    s = s / 100.0
    l = l / 100.0
    
    # Convert to RGB
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    
    # Convert to hex
    return "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))

# Custom OKLAB conversion functions
def linear_srgb_to_oklab(rgb):
    """Convert linear sRGB to OKLAB"""
    # Convert from sRGB to linear RGB
    linear_rgb = np.zeros_like(rgb)
    mask = rgb <= 0.04045
    linear_rgb[mask] = rgb[mask] / 12.92
    linear_rgb[~mask] = ((rgb[~mask] + 0.055) / 1.055) ** 2.4
    
    # Convert from linear RGB to OKLAB
    l = 0.4122214708 * linear_rgb[:, 0] + 0.5363325363 * linear_rgb[:, 1] + 0.0514459929 * linear_rgb[:, 2]
    m = 0.2119034982 * linear_rgb[:, 0] + 0.6806995451 * linear_rgb[:, 1] + 0.1073969566 * linear_rgb[:, 2]
    s = 0.0883024619 * linear_rgb[:, 0] + 0.2817188376 * linear_rgb[:, 1] + 0.6299787005 * linear_rgb[:, 2]
    
    l_ = np.cbrt(l)
    m_ = np.cbrt(m)
    s_ = np.cbrt(s)
    
    return np.column_stack([
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_
    ])

def oklab_to_linear_srgb(lab):
    """Convert OKLAB to linear sRGB"""
    l_ = lab[:, 0] + 0.3963377774 * lab[:, 1] + 0.2158037573 * lab[:, 2]
    m_ = lab[:, 0] - 0.1055613458 * lab[:, 1] - 0.0638541728 * lab[:, 2]
    s_ = lab[:, 0] - 0.0894841775 * lab[:, 1] - 1.2914855480 * lab[:, 2]
    
    l = l_ * l_ * l_
    m = m_ * m_ * m_
    s = s_ * s_ * s_
    
    linear_rgb = np.column_stack([
        +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
    ])
    
    # Clip to valid RGB range
    linear_rgb = np.clip(linear_rgb, 0, 1)
    
    # Convert from linear RGB to sRGB
    srgb = np.zeros_like(linear_rgb)
    mask = linear_rgb <= 0.0031308
    srgb[mask] = 12.92 * linear_rgb[mask]
    srgb[~mask] = 1.055 * (linear_rgb[~mask] ** (1/2.4)) - 0.055
    
    return srgb

# Color vision deficiency simulation
def cvd_convert(cols, severity=1, cvd_type="none"):
    """Convert colors to simulate color vision deficiency"""
    # Import needed libraries - included in function to avoid scoping issues
    from matplotlib.colors import to_rgb, to_hex
    import numpy as np
    from colorspacious import cspace_convert
    
    if cvd_type == "none":
        return cols
    
    # Convert hex to RGB
    rgb = np.array([to_rgb(c) for c in cols])
    
    # Map our cvd_type to the colorspacious naming
    cvd_type_map = {
        "deutan": "deuteranomaly",
        "protan": "protanomaly",
        "tritan": "tritanomaly"
    }
    
    # Use the correct CVD type
    if cvd_type in cvd_type_map:
        cvd_space = {"name": "sRGB1+CVD", 
                     "cvd_type": cvd_type_map[cvd_type], 
                     "severity": severity}
        
        try:
            # Convert to CVD space
            cvd_rgb = cspace_convert(rgb, "sRGB1", cvd_space)
            
            # Clip values to ensure they're in the valid range [0, 1]
            cvd_rgb = np.clip(cvd_rgb, 0, 1)
            
            # Convert back to hex
            return [to_hex(c) for c in cvd_rgb]
        except Exception as e:
            print(f"Warning: Error during color conversion ({e}). Returning original colors.")
            return cols
    else:
        return cols  # If unrecognized type, return original colors

# Flexible bounds function
def flex_bounds(x, xlim=(0, 1), rat=1.1):
    """Calculate flexible bounds based on data distribution"""
    range_x = np.max(x) - np.min(x)
    mean_x = np.mean(x)
    new_x = np.array([mean_x - range_x/2 * rat, mean_x + range_x/2 * rat])
    new_x[0] = max(xlim[0], new_x[0])
    new_x[1] = min(xlim[1], new_x[1])
    return new_x

# Color space ranges
cs_ranges = {
    "hsl": {
        "min": {"h": 0, "s": 0, "l": 0},
        "max": {"h": 360, "s": 100, "l": 100}
    },
    "lab": {
        "min": {"l": 0, "a": -128, "b": -128},
        "max": {"l": 100, "a": 127, "b": 127}
    },
    "lch": {
        "min": {"l": 0, "c": 0, "h": 0},
        "max": {"l": 100, "c": 140, "h": 360}
    },
    "oklab": {
        "min": {"l": 0, "a": -0.5, "b": -0.5},
        "max": {"l": 1, "a": 0.5, "b": 0.5}
    },
    "oklch": {
        "min": {"l": 0, "c": 0, "h": 0},
        "max": {"l": 1, "c": 0.4, "h": 360}
    }
}

# Rescale function
def rescale(df, color_space):
    """Rescale colors from [0,1] to their original ranges"""
    result = df.copy()
    for col in df.columns:
        min_val = cs_ranges[color_space]["min"][col]
        max_val = cs_ranges[color_space]["max"][col]
        result[col] = df[col] * (max_val - min_val) + min_val
    return result

def decode_color(hex_colors, to_space):
    """Convert hex colors to specified color space"""
    # Import locally to avoid scoping issues
    from matplotlib.colors import to_rgb
    import numpy as np
    import colorsys
    from colorspacious import cspace_convert
    
    # Convert hex to RGB first
    rgb = np.array([to_rgb(c) for c in hex_colors])
    
    if to_space == "hsl":
        # Convert RGB to HSL
        hsl_vals = np.array([colorsys.rgb_to_hls(r, g, b) for r, g, b in rgb])
        # Note: HLS to HSL reordering
        return pd.DataFrame({
            "h": hsl_vals[:, 0] * 360,
            "s": hsl_vals[:, 2] * 100,
            "l": hsl_vals[:, 1] * 100
        })
    elif to_space == "lab":
        # Convert RGB to LAB
        lab_vals = cspace_convert(rgb, "sRGB1", "CIELab")
        return pd.DataFrame({
            "l": lab_vals[:, 0],
            "a": lab_vals[:, 1],
            "b": lab_vals[:, 2]
        })
    elif to_space == "lch":
        # Convert RGB to LCH
        lab_vals = cspace_convert(rgb, "sRGB1", "CIELab")
        lch_vals = cspace_convert(lab_vals, "CIELab", "CIELCh")
        return pd.DataFrame({
            "l": lch_vals[:, 0],
            "c": lch_vals[:, 1],
            "h": lch_vals[:, 2]
        })
    elif to_space == "oklab":
        # Convert RGB to OKLAB using our custom function
        oklab_vals = linear_srgb_to_oklab(rgb)
        return pd.DataFrame({
            "l": oklab_vals[:, 0],
            "a": oklab_vals[:, 1],
            "b": oklab_vals[:, 2]
        })
    elif to_space == "oklch":
        # Convert RGB to OKLCH (via OKLAB)
        oklab_vals = linear_srgb_to_oklab(rgb)
        # Convert rectangular to polar coordinates for a/b channels
        l_vals = oklab_vals[:, 0]
        c_vals = np.sqrt(oklab_vals[:, 1]**2 + oklab_vals[:, 2]**2)
        h_vals = np.arctan2(oklab_vals[:, 2], oklab_vals[:, 1]) * 180 / np.pi
        # Ensure h is in range [0, 360]
        h_vals = np.mod(h_vals + 360, 360)
        return pd.DataFrame({
            "l": l_vals,
            "c": c_vals,
            "h": h_vals
        })
    else:
        raise ValueError(f"Unsupported color space: {to_space}")

def encode_color(df, from_space):
    """Convert colors from specified color space to hex"""
    # Import locally to avoid scoping issues
    from matplotlib.colors import to_hex
    import numpy as np
    import colorsys
    from colorspacious import cspace_convert
    
    values = df.values if isinstance(df, pd.DataFrame) else df
    
    if from_space == "hsl":
        # HSL to RGB
        rgb_vals = np.array([
            colorsys.hls_to_rgb(h/360, l/100, s/100)
            for h, s, l in values[:, [0, 1, 2]]
        ])
    elif from_space == "lab":
        # LAB to RGB
        rgb_vals = cspace_convert(values, "CIELab", "sRGB1")
    elif from_space == "lch":
        # LCH to RGB (via LAB)
        lab_vals = cspace_convert(values, "CIELCh", "CIELab")
        rgb_vals = cspace_convert(lab_vals, "CIELab", "sRGB1")
    elif from_space == "oklab":
        # OKLAB to RGB using our custom function
        rgb_vals = oklab_to_linear_srgb(values)
    elif from_space == "oklch":
        # OKLCH to RGB (via OKLAB)
        l_vals = values[:, 0]
        c_vals = values[:, 1]
        h_vals = values[:, 2] * np.pi / 180  # convert to radians
        
        # Convert polar to rectangular
        a_vals = c_vals * np.cos(h_vals)
        b_vals = c_vals * np.sin(h_vals)
        
        oklab_vals = np.column_stack([l_vals, a_vals, b_vals])
        rgb_vals = oklab_to_linear_srgb(oklab_vals)
    else:
        raise ValueError(f"Unsupported color space: {from_space}")
    
    # Clip RGB values to valid range [0, 1]
    rgb_vals = np.clip(rgb_vals, 0, 1)
    
    # Convert to hex
    return [to_hex(rgb) for rgb in rgb_vals]

# Implementation of CIEDE2000 color difference formula
def deltaE_ciede2000(lab1, lab2, kL=1, kC=1, kH=1):
    """
    Calculate the CIEDE2000 color difference between two CIELAB colors.
    """
    # Extract coordinates
    L1, a1, b1 = lab1
    L2, a2, b2 = lab2
    
    # Calculate C1, C2 (chroma)
    C1 = np.sqrt(a1**2 + b1**2)
    C2 = np.sqrt(a2**2 + b2**2)
    
    # Calculate mean C
    Cbar = (C1 + C2) / 2
    
    # Calculate G
    G = 0.5 * (1 - np.sqrt(Cbar**7 / (Cbar**7 + 25**7)))
    
    # Calculate a'
    a1p = (1 + G) * a1
    a2p = (1 + G) * a2
    
    # Calculate C'
    C1p = np.sqrt(a1p**2 + b1**2)
    C2p = np.sqrt(a2p**2 + b2**2)
    
    # Calculate h' (hue angle)
    h1p = np.arctan2(b1, a1p) % (2 * np.pi)
    h2p = np.arctan2(b2, a2p) % (2 * np.pi)
    
    # Calculate ΔL', ΔC', ΔH'
    dLp = L2 - L1
    dCp = C2p - C1p
    
    # Calculate ΔH'
    dhp = h2p - h1p
    if dhp > np.pi:
        dhp -= 2 * np.pi
    elif dhp < -np.pi:
        dhp += 2 * np.pi
    
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(dhp / 2)
    
    # Calculate mean L', C', h'
    Lbar = (L1 + L2) / 2
    Cbar = (C1p + C2p) / 2
    
    # Calculate mean h'
    if C1p * C2p == 0:
        hbar = h1p + h2p
    else:
        if abs(h1p - h2p) <= np.pi:
            hbar = (h1p + h2p) / 2
        else:
            if h1p + h2p < 2 * np.pi:
                hbar = (h1p + h2p + 2 * np.pi) / 2
            else:
                hbar = (h1p + h2p - 2 * np.pi) / 2
    
    # Calculate T
    T = 1 - 0.17 * np.cos(hbar - np.pi/6) + 0.24 * np.cos(2*hbar) + \
        0.32 * np.cos(3*hbar + np.pi/30) - 0.2 * np.cos(4*hbar - 63*np.pi/180)
    
    # Calculate ΔΘ
    dtheta = 30 * np.exp(-((hbar*180/np.pi - 275) / 25)**2)
    
    # Calculate RC
    RC = 2 * np.sqrt(Cbar**7 / (Cbar**7 + 25**7))
    
    # Calculate SL, SC, SH
    SL = 1 + (0.015 * (Lbar - 50)**2) / np.sqrt(20 + (Lbar - 50)**2)
    SC = 1 + 0.045 * Cbar
    SH = 1 + 0.015 * Cbar * T
    
    # Calculate RT
    RT = -np.sin(2 * dtheta * np.pi / 180) * RC
    
    # Calculate color difference
    dE00 = np.sqrt(
        (dLp / (kL * SL))**2 +
        (dCp / (kC * SC))**2 +
        (dHp / (kH * SH))**2 +
        RT * (dCp / (kC * SC)) * (dHp / (kH * SH))
    )
    
    return dE00

# Color comparison function
def compare_color(from_colors, to_colors=None, method="cie2000", from_space="lab"):
    """Calculate color difference between colors"""
    # Import locally to avoid scoping issues
    from matplotlib.colors import to_rgb
    import numpy as np
    from colorspacious import cspace_convert
    
    # Convert to lab space for comparison if needed
    if from_space != "lab":
        # First convert to sRGB
        if isinstance(from_colors, pd.DataFrame):
            from_colors_rescaled = rescale(from_colors, from_space)
            from_rgb = np.array([to_rgb(c) for c in encode_color(from_colors_rescaled, from_space)])
        else:
            from_rgb = np.array([to_rgb(c) for c in from_colors])
        
        # Then to LAB
        from_lab = cspace_convert(from_rgb, "sRGB1", "CIELab")
        
        if to_colors is not None:
            if isinstance(to_colors, pd.DataFrame):
                to_colors_rescaled = rescale(to_colors, from_space)
                to_rgb = np.array([to_rgb(c) for c in encode_color(to_colors_rescaled, from_space)])
            else:
                to_rgb = np.array([to_rgb(c) for c in to_colors])
            to_lab = cspace_convert(to_rgb, "sRGB1", "CIELab")
    else:
        from_lab = from_colors.values if isinstance(from_colors, pd.DataFrame) else from_colors
        if to_colors is not None:
            to_lab = to_colors.values if isinstance(to_colors, pd.DataFrame) else to_colors
    
    # Calculate distances
    if to_colors is None:
        # Pairwise distances within from_colors
        n = len(from_lab)
        dists = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                if method == "cie76":
                    # Euclidean distance
                    dist = np.sqrt(np.sum((from_lab[i] - from_lab[j])**2))
                elif method == "cie2000":
                    # Delta E 2000
                    dist = deltaE_ciede2000(from_lab[i], from_lab[j])
                else:
                    raise ValueError(f"Unsupported distance method: {method}")
                
                dists[i, j] = dists[j, i] = dist
        
        return dists
    else:
        # Distances between from_colors and to_colors
        nf = len(from_lab)
        nt = len(to_lab)
        dists = np.zeros((nf, nt))
        
        for i in range(nf):
            for j in range(nt):
                if method == "cie76":
                    # Euclidean distance
                    dists[i, j] = np.sqrt(np.sum((from_lab[i] - to_lab[j])**2))
                elif method == "cie2000":
                    # Delta E 2000
                    dists[i, j] = deltaE_ciede2000(from_lab[i], to_lab[j])
                else:
                    raise ValueError(f"Unsupported distance method: {method}")
        
        return dists

# Objective function for optimization
def mean_dist_func(par, curr_cols, bounds_sc, bounds_l, color_space, 
                  dist_alg="cie2000", colorblind_safe=True, 
                  colorblind_weights={"none": 6, "deutan": 6, "protan": 2, "tritan": 0.1},
                  init_weights=None, init_colors=None, return_info=False):
    """Objective function to maximize harmonic mean of pairwise color distances"""
    cn = curr_cols.columns.tolist()
    
    # Transform parameters to [0,1]
    m = par.reshape(-1, len(cn))
    nncols = m.shape[0]
    m = pd.DataFrame(m, columns=cn)
    
    # For luminosity
    if "l" in cn:
        if nncols > 1:
            m.loc[1:, "l"] = np.exp(m.loc[1:, "l"].values)
            m["l"] = m["l"].cumsum()
        m["l"] = invlogit(m["l"].values) * (bounds_l[1] - bounds_l[0]) + bounds_l[0]
    
    # For chroma / saturation
    sc_cols = [col for col in cn if col in ["s", "c"]]
    if sc_cols:
        for col in sc_cols:
            m[col] = invlogit(m[col].values) * (bounds_sc[1] - bounds_sc[0]) + bounds_sc[0]
    
    # For remaining channels
    other_channels = [col for col in cn if col not in ["s", "c", "l"]]
    if other_channels:
        for col in other_channels:
            m[col] = invlogit(m[col].values)
    
    # Convert back to real units and recover hex codes
    new_hex = encode_color(rescale(m, color_space), color_space)
    curr_hex = encode_color(rescale(curr_cols, color_space), color_space)
    
    # Convert to colorblind-space and compute pairwise distances
    cvd_states = ["deutan", "protan", "tritan", "none"] if colorblind_safe else ["none"]
    
    d = {cvd_state: 0 for cvd_state in cvd_states}
    for cvd_state in cvd_states:
        # Simulate colorblindness
        nh = cvd_convert(new_hex, cvd_type=cvd_state)
        ch = cvd_convert(curr_hex, cvd_type=cvd_state)
        
        # Get back into appropriate colorspace
        cvals = decode_color(ch, color_space)
        nvals = decode_color(nh, color_space)
        
        # Find distances between colors
        # Old colors to new colors
        dists1 = compare_color(cvals, nvals, method=dist_alg, from_space=color_space)
        dists1 = dists1.flatten()  # Flatten matrix to array
        
        # New colors to new colors
        if len(nh) > 1:
            dists2 = compare_color(nvals, method=dist_alg, from_space=color_space)
            # Get upper triangle
            dists2 = dists2[np.triu_indices(len(nh), k=1)]
            dists = np.concatenate([dists1, dists2])
        else:
            dists = dists1
        
        # Compute harmonic mean of distances
        hmean_dist = len(dists) / np.sum(1.0 / np.maximum(dists, 1e-10))
        d[cvd_state] = hmean_dist
    
    # Compute weighted harmonic mean
    wd = sum(d[cvd_state] * colorblind_weights.get(cvd_state, 0) for cvd_state in cvd_states)
    
    # Also compute a penalty term so values can't rocket away
    penalty = np.sum(par**2)
    
    # Are we trying to also stick to some initial colors?
    if init_weights is not None and init_colors is not None:
        init_cvals = decode_color(init_colors, color_space)
        dists_to_init = np.array([
            compare_color(init_cvals.iloc[i:i+1], nvals.iloc[i:i+1], method=dist_alg, from_space=color_space)[0][0]
            for i in range(nncols)
        ])
        weighted_dists_to_init = np.sum(dists_to_init * init_weights)
        penalty += weighted_dists_to_init
    
    # Do we want to return the objective, or other relevant metadata?
    if return_info:
        out = {
            "wd": wd,
            "new_hex": new_hex
        }
    else:
        # Since optimize minimizes but we want to maximize the distance, multiply by -1
        out = -wd + penalty
    
    return out

# Optimization function - SERIAL VERSION
def optimize_colors(
    curr_cols, 
    n_cols_to_add=1,
    n_optim_runs=20,
    color_space="oklab",
    bounds_sc=(0, 1),
    bounds_l=(0, 1),
    colorblind_safe=True,
    colorblind_weights={"none": 6, "deutan": 6, "protan": 2, "tritan": 0.1}
):
    """Run optimization to find optimal new colors"""
    np.random.seed(42)  # For reproducibility
    
    # Initialize parameter vector for each optimization run
    starts = [np.random.normal(0, 1, n_cols_to_add * len(curr_cols.columns)) for _ in range(n_optim_runs)]
    
    # Run optimization serially
    results = []
    for i, start_params in enumerate(starts):
        print(f"Running optimization {i+1}/{n_optim_runs}...")
        try:
            result = optimize.minimize(
                mean_dist_func,
                x0=start_params,
                args=(curr_cols, bounds_sc, bounds_l, color_space, "cie2000", colorblind_safe, colorblind_weights),
                method='BFGS'
            )
            results.append((result.fun, result.x))
        except Exception as e:
            print(f"Optimization run {i+1} failed: {e}")
            # Add a dummy result with high objective value
            results.append((1e10, np.zeros_like(start_params)))
    
    if not results:
        raise ValueError("All optimization runs failed")
    
    # Find optimal result
    best_idx = np.argmin([r[0] for r in results])
    raw_pars = results[best_idx][1]
    
    # Generate new colors
    new_colors = mean_dist_func(
        par=raw_pars,
        curr_cols=curr_cols,
        bounds_sc=bounds_sc,
        bounds_l=bounds_l,
        color_space=color_space,
        colorblind_weights=colorblind_weights,
        colorblind_safe=colorblind_safe,
        return_info=True
    )
    
    return new_colors["new_hex"]

# Visualization functions
def plot_color_addition(cols, new_cols, plot_cbs=["none", "deutan", "protan", "tritan"], 
                        colorwheel_space="hsl", file_path="color_addition.png"):
    """
    Plot current and new colors, with simulations for different types of color blindness
    """
    # Import locally to avoid scoping issues
    from matplotlib.colors import to_rgb
    import numpy as np
    from colorspacious import cspace_convert
    
    n_cbs = len(plot_cbs)
    fig = plt.figure(figsize=(7.5 * n_cbs, 10), dpi=100)
    gs = gridspec.GridSpec(2, n_cbs)
    
    for cb_idx, plot_cb in enumerate(plot_cbs):
        # Convert colors for colorblindness simulation
        cb_cols = cvd_convert(cols, cvd_type=plot_cb)
        cb_new_cols = cvd_convert(new_cols, cvd_type=plot_cb)
        
        # Plot color bars
        ax1 = fig.add_subplot(gs[0, cb_idx])
        ax1.set_xlim(0, 5)
        ax1.set_ylim(0, 1)
        ax1.axis('off')
        
        max_col_num = max(len(cols), len(new_cols))
        col_height = 0.9 / max_col_num
        
        # Plot current colors
        for i, col in enumerate(cb_cols):
            rect = patches.Rectangle(
                (0.5, sum([col_height] * i)), 
                1, col_height, 
                color=col
            )
            ax1.add_patch(rect)
            # Determine text color (black or white)
            lab = cspace_convert(to_rgb(col), "sRGB1", "CIELab")
            text_color = 'white' if lab[0] < 50 else 'black'
            ax1.text(1, sum([col_height] * i) + col_height/2, cols[i], 
                     color=text_color, ha='center', va='center')
        
        # Plot new colors
        for i, col in enumerate(cb_new_cols):
            rect = patches.Rectangle(
                (2, sum([col_height] * i)), 
                1, col_height, 
                color=col
            )
            ax1.add_patch(rect)
            # Determine text color
            lab = cspace_convert(to_rgb(col), "sRGB1", "CIELab")
            text_color = 'white' if lab[0] < 50 else 'black'
            ax1.text(2.5, sum([col_height] * i) + col_height/2, new_cols[i], 
                     color=text_color, ha='center', va='center')
        
        # Add labels
        ax1.text(1, len(cols) * col_height, "current colors = O", fontsize=15, ha='center', va='bottom')
        ax1.text(2.5, len(new_cols) * col_height, "new colors = □", fontsize=15, ha='center', va='bottom')
        
        # Add luminosity scale - FIXED to use hex instead of HSL strings
        nsegs = 50
        for i in range(nsegs):
            lum = i / nsegs * 100
            # Use our helper function to convert HSL to hex that matplotlib can understand
            color = hsl_to_hex(0, 0, lum)
            rect = patches.Rectangle(
                (3.9, i/nsegs * 0.5), 
                0.2, 0.5/nsegs, 
                color=color
            )
            ax1.add_patch(rect)
        
        # Add border for luminosity scale
        rect = patches.Rectangle(
            (3.9, 0), 
            0.2, 0.5 + 1/nsegs, 
            fill=False,
            edgecolor='black'
        )
        ax1.add_patch(rect)
        
        # Add scale dots
        for i in range(10):
            lum_val = i / 9 * 0.5
            ax1.scatter(3.8 - i/90, lum_val, marker='o', 
                      s=(i/9*100/50 + 0.5)*100, edgecolor='black', facecolor='none')
            ax1.scatter(4.2 + i/90, lum_val, marker='s', 
                      s=(i/9*100/50 + 0.5)*100, edgecolor='black', facecolor='none')
        
        ax1.text(4, 0.5 + 1/nsegs, "luminosity", fontsize=15, ha='center', va='bottom')
        
        # Plot color wheel
        ax2 = fig.add_subplot(gs[1, cb_idx])
        ax2.set_xlim(-1, 1)
        ax2.set_ylim(-1, 1)
        ax2.axis('off')
        ax2.set_aspect('equal')
        
        # Create color wheel
        n_slices_t = 60
        n_slices_r = 20
        degrees = np.linspace(0, 2*np.pi, n_slices_t+1)
        dd = np.diff(degrees)[0]/2
        radii = np.linspace(0, 1, n_slices_r+1)
        rd = np.diff(radii)[0]/2
        
        for ts in range(n_slices_t):
            for rs in range(n_slices_r):
                # Get coordinates for the patch
                coords = np.array([
                    polar2cart(degrees[ts], radii[rs]),
                    polar2cart(degrees[ts], radii[rs+1]),
                    polar2cart(degrees[ts+1], radii[rs+1]),
                    polar2cart(degrees[ts+1], radii[rs])
                ])
                
                # Get center hue/chroma of this slice
                hue_deg = (degrees[ts] + dd) / (2 * np.pi) * 360
                chroma = (radii[rs] + rd)
                
                # Create color for wheel patch
                fixed_lightness = 0.5 if colorwheel_space == "hsl" else 0.75
                
                if colorwheel_space == "hsl":
                    wheel_col = encode_color(
                        pd.DataFrame({
                            "h": [hue_deg],
                            "s": [chroma * cs_ranges["hsl"]["max"]["s"]],
                            "l": [fixed_lightness * cs_ranges["hsl"]["max"]["l"]]
                        }), "hsl")[0]
                elif colorwheel_space == "lch":
                    wheel_col = encode_color(
                        pd.DataFrame({
                            "l": [fixed_lightness * cs_ranges["lch"]["max"]["l"]],
                            "c": [chroma * cs_ranges["lch"]["max"]["c"]],
                            "h": [hue_deg]
                        }), "lch")[0]
                elif colorwheel_space == "oklch":
                    wheel_col = encode_color(
                        pd.DataFrame({
                            "l": [fixed_lightness],
                            "c": [chroma],
                            "h": [hue_deg]
                        }), "oklch")[0]
                
                # Convert for colorblindness
                wheel_col = cvd_convert([wheel_col], cvd_type=plot_cb)[0]
                
                # Add patch to color wheel
                poly = patches.Polygon(coords, closed=True, fill=True, 
                                      edgecolor=None, facecolor=wheel_col)
                ax2.add_patch(poly)
        
        # Decode colors to chosen colorwheel space
        colorwheel_cols = decode_color(cols, colorwheel_space)
        colorwheel_new_cols = decode_color(new_cols, colorwheel_space)
        
        # Get coordinate mapping for current and new colors
        s_or_c = "s" if "s" in colorwheel_cols.columns else "c"
        
        # Current colors coordinates
        coords_curr_cols = np.array([
            polar2cart(row["h"] / 360 * 2 * np.pi, 
                      row[s_or_c] / cs_ranges[colorwheel_space]["max"][s_or_c])
            for _, row in colorwheel_cols.iterrows()
        ])
        
        # New colors coordinates
        coords_new_cols = np.array([
            polar2cart(row["h"] / 360 * 2 * np.pi, 
                      row[s_or_c] / cs_ranges[colorwheel_space]["max"][s_or_c])
            for _, row in colorwheel_new_cols.iterrows()
        ])
        
        # Choose point sizes based on lightness
        l_min = cs_ranges[colorwheel_space]["min"]["l"]
        l_max = cs_ranges[colorwheel_space]["max"]["l"]
        
        l_old_norm = (colorwheel_cols["l"] - l_min) / (l_max - l_min)
        l_new_norm = (colorwheel_new_cols["l"] - l_min) / (l_max - l_min)
        
        cex_old = 0.5 + 2 * l_old_norm
        cex_new = 0.5 + 2 * l_new_norm
        
        # Plot original and new colors
        for i, (coord, col, size) in enumerate(zip(coords_curr_cols, cvd_convert(cols, cvd_type=plot_cb), cex_old)):
            ax2.scatter(coord[0], coord[1], s=size*100, marker='o', 
                       edgecolor='black', facecolor=col)
        
        for i, (coord, col, size) in enumerate(zip(coords_new_cols, cvd_convert(new_cols, cvd_type=plot_cb), cex_new)):
            ax2.scatter(coord[0], coord[1], s=size*100, marker='s', 
                       edgecolor='black', facecolor=col)
        
        # Add legend
        legend_text = f"colorwheel in\n{colorwheel_space} colorspace"
        if plot_cb != "none":
            legend_text += f"\n({plot_cb}-type)"
        ax2.text(-0.9, 0.9, legend_text, fontsize=12, va='top')
    
    plt.tight_layout()
    plt.savefig(file_path)
    plt.close()
    return file_path

# Main function to generate colors
def generate_closest_colors(
    colors=None,
    color_space="oklab",
    colorblind_safe=True,
    colorblind_weights={"none": 6, "deutan": 6, "protan": 2, "tritan": 0.1},
    brewer_cols=True,
    rainbow_circle_cols=False,
    constrain_to_similar_colors=True,
    fb_rat=0.9,
    n_rainbow=6,
    n_optim_runs=100,
    n_cols_to_add=1,
    colorwheel_space="hsl",
    plot_cbs=["none", "deutan", "protan", "tritan"],
    output_path="color_addition.png"
):
    """
    Main function to generate optimal color palettes
    
    Parameters:
    -----------
    colors : list, optional
        List of hex color codes to start with
    color_space : str, default="oklab"
        Color space to optimize in, one of "hsl", "lab", "lch", "oklab", "oklch"
    colorblind_safe : bool, default=True
        Whether to account for colorblindness
    colorblind_weights : dict
        Weights for different types of colorblindness
    brewer_cols : bool, default=True
        Whether to use ColorBrewer palette as starting point
    rainbow_circle_cols : bool, default=False
        Whether to use evenly spaced hues as starting point
    constrain_to_similar_colors : bool, default=True
        Whether to constrain new colors to similar ranges as existing ones
    fb_rat : float, default=0.9
        How much wider the constraint should be than the quantiles
    n_rainbow : int, default=6
        Number of colors if using rainbow_circle_cols
    n_optim_runs : int, default=100
        Number of optimization runs
    n_cols_to_add : int, default=1
        Number of colors to add
    colorwheel_space : str, default="hsl"
        Color space for visualization
    plot_cbs : list
        Types of colorblindness to simulate in visualization
    output_path : str, default="color_addition.png"
        Path to save output image
        
    Returns:
    --------
    dict
        Contains 'original_colors', 'new_colors', and 'plot_path'
    """
    # Import locally to avoid scoping issues
    from matplotlib.colors import to_hex
    import numpy as np
    
    # Generate initial colors
    if colors is not None:
        cols = colors
    elif brewer_cols:
        # Using matplotlib's Set1 which is similar to RColorBrewer::brewer.pal(8, "Set1")
        cols = plt.cm.Set1.colors[:8]
        cols = [to_hex(c) for c in cols]
    elif rainbow_circle_cols:
        # Create evenly spaced hues
        hues = np.linspace(0, 360, n_rainbow + 1)[:-1]
        cols = [hsl_to_hex(h, 100, 50) for h in hues]  # Use our helper function instead
    else:
        raise ValueError("Must provide colors or set brewer_cols or rainbow_circle_cols to True")
    
    # Decode colors to specified colorspace
    dat = decode_color(cols, color_space)
    
    # Calculate bounds if necessary
    if constrain_to_similar_colors:
        # For the saturation or chroma channel, if one exists
        s_or_c_col = next((col for col in dat.columns if col in ["s", "c"]), None)
        if s_or_c_col:
            bounds_sc = flex_bounds(
                np.quantile(dat[s_or_c_col], [0.05, 0.95]), 
                rat=fb_rat
            )
        else:
            bounds_sc = np.array([0, 1])
        
        # For the luminosity channel
        bounds_l = flex_bounds(
            np.quantile(dat["l"], [0.1, 0.9]), 
            rat=fb_rat
        )
    else:
        bounds_sc = np.array([0, 1])
        bounds_l = np.array([0, 1])
    
    # Optimize colors
    new_cols = optimize_colors(
        curr_cols=dat,
        n_cols_to_add=n_cols_to_add,
        n_optim_runs=n_optim_runs,
        color_space=color_space,
        bounds_sc=bounds_sc,
        bounds_l=bounds_l,
        colorblind_safe=colorblind_safe,
        colorblind_weights=colorblind_weights
    )
    
    # Plot results
    plot_path = plot_color_addition(
        cols=cols,
        new_cols=new_cols,
        plot_cbs=plot_cbs,
        colorwheel_space=colorwheel_space,
        file_path=output_path
    )
    
    return {
        "original_colors": cols,
        "new_colors": new_cols,
        "plot_path": plot_path
    }

# Example usage
if __name__ == "__main__":
    # Set parameters (matching R script defaults)
    params = {
        "color_space": "oklab",  # Using oklab as in R script
        "colorblind_safe": True,
        "colorblind_weights": {"none": 6, "deutan": 6, "protan": 2, "tritan": 0.1},
        "brewer_cols": True,
        "rainbow_circle_cols": False,
        "constrain_to_similar_colors": True,
        "fb_rat": 0.9,
        "n_rainbow": 6,
        "n_optim_runs": 1,  # Reduced for quicker results
        "n_cols_to_add": 3,
        "colorwheel_space": "hsl",
        "plot_cbs": ["none", "deutan", "protan", "tritan"],
        "output_path": "color_addition.png"
    }
    
    # Example with custom colors
    custom_colors = {
        "Transcriptomics": "#4477AA",
        "Metabolomics": "#6D4B08",
        "Proteomics": "#228833",
        "Phosphoproteomics": "#F3A02B",
        "ATAC": "#882255",
        "Methylation": "#D687B5"
    }
    
    result = generate_closest_colors(
        colors=list(custom_colors.values()),
        **params
    )
    
    print("Original colors:", result["original_colors"])
    print("New colors:", result["new_colors"])
    print("Plot saved to:", result["plot_path"])