"""Unified builder for the 9-block blade-section grid topology.

Reproduces the OGV55 midspan working baseline (1 O-grid + 3-block top passage +
3-block bot passage + 2 LE/TE wrap blocks; ~16500 grid points, 0 folds across
all blocks). Consolidates these previously separate scripts:

    build_ogrid_3.py
    build_top_passage_3.py
    build_bot_passage_3.py
    build_middle_passage_3.py

Input:  one T-Blade3 blade-section file (the "blade.X.X.SECTION" output from
        T-Blade3, two-column m'/theta with a 2-line header — line 1 is the case
        name, line 2 has chord/pitch in positions 3 and 4).
Output: <output_dir>/{ogrid_3, top_passage_3/{left,center,right},
        bot_passage_3/{left,center,right}, middle_passage_3/{left,right}},
        each containing x_boundaries.dat, y_boundaries.dat, case.in — ready
        to run with `gridexec case.in`.

Tunables are exposed at the top of the file. Defaults reproduce the OGV55
midspan working values (lambda_max per block, INFLATE_DIST, M_LEFT/RIGHT_END,
DS_OGRID_FIRST, anchor offsets). For a different blade section the most
likely parameters needing tuning are M_LEFT_END, M_RIGHT_END, INFLATE_DIST_*,
and the per-block lambda_max values (see project_blade_section_lambda_sweep.md).

Usage:
    python build_9block_topology.py BLADE_PATH OUTPUT_DIR
"""

import argparse
import os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


# ============================================================================
# TUNABLES (defaults reproduce OGV55 midspan working baseline)
# ============================================================================

# --- Geometry ---
OFFSET_FRAC          = 0.02      # ogrid contour offset = OFFSET_FRAC * chord
INFLATE_DIST_TOP     = 0.10      # offset of inflated upperSS lid
INFLATE_DIST_BOT     = 0.10      # offset of inflated lowerPS lid (must stay below R_min ~ 0.131 for OGV55)
M_LEFT_END           = -0.35     # left wall in m' (vertical wall at xi=0 of LEFT blocks)
M_RIGHT_END          =  0.46     # right wall in m'

# --- Anchor selection (PER-SIDE ADAPTIVE) ---
# Default symmetric anchors are tried first (A_uLE = N/2 - LE_OFFSET_DEFAULT,
# A_lLE = N/2 + LE_OFFSET_DEFAULT). The LE arc theta range is checked against
# the anchor thetas; if either side's margin is more negative than
# LE_INTERSECTION_TOL_FRAC * chord, ONLY that side is adapted by walking to
# the nearest local theta extremum. This preserves the OGV55 midspan baseline
# byte-identically (its tiny -4e-5 dip is below tolerance) while fixing
# heavily cambered sections like hub. After adaptation, the snapped anchor is
# nudged to keep LE arc length odd (so TE_OFFSET = (arc_len-1)/2 is integer).
LE_OFFSET_DEFAULT       = 20      # default anchor offset from i_LE_topo
LE_OFFSET_MIN           = 5       # adapt search starts here
LE_OFFSET_MAX           = 40      # adapt search bound
LE_INTERSECTION_TOL_FRAC = 1e-3   # tolerance as fraction of chord
LE_GEOM_TOPO_TOL        = 10      # warn if |argmin(o_m) - N/2| exceeds this

# Sharp ds-ratio jumps in the offset-contour near the TE-tip (T-Blade3 puts
# fine clustering at the TE rounded tip, then transitions to coarser SS/PS
# nodes). When A_uTE / A_lTE land at the transition node, inflate_arc with
# INFLATE_DIST=0.10 amplifies the local metric jump and folds the centers.
# Walk te_offset outward (parity-preserving, +2 per step) to skip the
# transition; cap at TE_OFFSET_WALK_MAX. The walker requires not just A_uTE
# itself but also TE_KINK_BUFFER consecutive nodes interior to it (toward the
# closure) to be smooth — this puts the kink mid-wall in middle.right rather
# than adjacent to its corner, where corner clustering can't absorb it.
TE_DS_RATIO_TOL         = 1.4     # max(ds_left, ds_right)/min < this is "smooth"
TE_OFFSET_WALK_MAX      = 15      # walk up to this many nodes past LE-derived te_offset
TE_KINK_BUFFER          = 2       # require this many smooth nodes interior to A_uTE

# --- Block grid sizes ---
N_OGRID_RADIAL       = 5         # eta points in O-grid block (xi count = blade pts after closure)
DS_OGRID_FIRST       = 0.0018    # first radial spacing in O-grid (eta=0 -> eta=1)
N_X_LR               = 31        # xi-pts in LEFT/RIGHT L-blocks of top/bot passages
                                 # (n_y_passage is now derived from the adaptive
                                 # LE arc length; flows through find_anchors.)

# --- Shared L-block clustering (LE/TE corner first-segment scales) ---
# Working OGV55 values. Independently re-tune per blade section if needed
# (smaller -> tighter clustering near the LE/TE corner). The build scripts'
# original comments described these as "avg of upperSS + lowerPS smoothed_ds
# at the LE corner" / "TE-side smoothed_ds (top only, bot is cusp artifact)",
# but the recorded values do not reproduce from a fresh smoothed_ds(k=5)
# evaluation against the current ogrid contour — they were hand-tuned. The
# OGV55 baseline (0 folds, 9 blocks) was achieved with these.
SHARED_LE_DS         = 0.002635  # cluster scale at LE corner of L-blocks
SHARED_TE_DS         = 0.000716  # cluster scale at TE corner of L-blocks

# --- Per-block Newton tunables (working OGV55 values) ---
LAMBDA_TOP   = {"left": 0.10, "center": 0.05, "right": 0.10}
LAMBDA_BOT   = {"left": 0.10, "center": 0.02, "right": 0.05}
LAMBDA_MID   = {"left": 0.10, "right": 0.10}
LAMBDA_OGRID = 0.10              # used by the psn_lim solver in the O-grid block

RAMP_ITERS_NEWTON = 20
NITER_NEWTON      = 200
RAMP_ITERS_OGRID  = 100
INNER_NITER_OGRID = 10
NITER_OGRID       = 500

TOL_STEP  = 1.0e-8
TOL_RESID = 1.0e-6


# ============================================================================
# T-Blade3 input
# ============================================================================

def read_tblade3_blade(path):
    """Read a T-Blade3 blade.X.X.SECTION file.

    Returns:
        m, t   : (N,) arrays of unique blade points (closing duplicate dropped)
        chord  : float (header field 3)
        pitch  : float (header field 4)
    """
    with open(path) as f:
        f.readline()                                          # case name
        header = [float(x) for x in f.readline().split()]
    chord = header[3]
    pitch = header[4]
    raw = np.loadtxt(path, skiprows=2)
    if raw.shape[1] != 2:
        raise ValueError(f"expected 2 columns, got {raw.shape[1]}")
    m = raw[:-1, 0]                                           # drop closing duplicate
    t = raw[:-1, 1]
    return m, t, chord, pitch


# ============================================================================
# Helpers
# ============================================================================

def signed_area(closed):
    """Signed polygon area; > 0 for CCW orientation."""
    return 0.5 * np.sum(closed[:-1, 0] * closed[1:, 1] -
                        closed[1:, 0] * closed[:-1, 1])


def build_ogrid_contour(m, t, offset_dist):
    """Offset the closed blade outward by offset_dist along the local outward
    normal (computed via index-parametrised periodic cubic spline). Returns
    the (N+1, 2) closed contour."""
    N = len(m)
    idx = np.arange(N + 1)
    cs_m = CubicSpline(idx, np.r_[m, m[0]], bc_type='periodic')
    cs_t = CubicSpline(idx, np.r_[t, t[0]], bc_type='periodic')
    dm = cs_m(np.arange(N), 1)
    dt = cs_t(np.arange(N), 1)
    mag = np.hypot(dm, dt)
    tx, ty = dm / mag, dt / mag

    # Pick the normal pointing OUTWARD (away from the centroid).
    mc, tc = m.mean(), t.mean()
    n1x, n1y = ty, -tx
    if (n1x * (m - mc) + n1y * (t - tc)).mean() > 0:
        nx_n, ny_n = n1x, n1y
    else:
        nx_n, ny_n = -ty, tx
    o_m = m + offset_dist * nx_n
    o_t = t + offset_dist * ny_n
    return np.column_stack([np.r_[o_m, o_m[0]], np.r_[o_t, o_t[0]]])


def find_anchors(ogrid_closed, chord, le_offset_default, le_offset_min,
                 le_offset_max, le_intersection_tol_frac, le_geom_topo_tol):
    """Per-side adaptive anchor selection.

    Returns (i_LE_geom, i_LE_topo, A_uTE, A_uLE, A_lLE, A_lTE, N, n_y_passage).

    Starts with default symmetric anchors (offset = le_offset_default). Checks
    whether the LE arc theta range fits inside [theta(A_lLE), theta(A_uLE)].
    For each side that violates by more than chord*tol_frac, walks the
    contour to the nearest local theta extremum and snaps to the anchor with
    parity preserving an odd LE arc length. Sides that pass the check stay
    untouched — keeps the OGV55 midspan baseline byte-identical."""
    if (ogrid_closed.shape[0] - 1) % 2 != 0:
        raise RuntimeError("expected even number of unique contour points "
                           "(T-Blade3 convention).")
    N = ogrid_closed.shape[0] - 1
    i_LE_topo = N // 2
    o_m = ogrid_closed[:N, 0]
    o_t = ogrid_closed[:N, 1]
    i_LE_geom = int(np.argmin(o_m))
    if abs(i_LE_geom - i_LE_topo) > le_geom_topo_tol:
        print(f"  WARNING: argmin(o_m)={i_LE_geom} is more than "
              f"{le_geom_topo_tol} nodes from N/2={i_LE_topo} — the blade "
              f"distribution may not follow the T-Blade3 N/2-per-side "
              f"convention; LE arc may not contain the geometric LE.")

    A_uLE = i_LE_topo - le_offset_default
    A_lLE = i_LE_topo + le_offset_default
    tol = le_intersection_tol_frac * chord
    arc = o_t[A_uLE:A_lLE + 1]
    margin_top = o_t[A_uLE] - arc.max()
    margin_bot = arc.min() - o_t[A_lLE]

    adapted = []

    # Bot adapt: walk forward from i_LE_topo + le_offset_min, find local min
    if margin_bot < -tol:
        k = i_LE_topo + le_offset_min
        while k < i_LE_topo + le_offset_max:
            if o_t[(k + 1) % N] > o_t[k]:
                break
            k += 1
        else:
            print(f"  WARNING: no local PS theta min within {le_offset_max} "
                  f"nodes; clamping at le_offset_max.")
        # Snap to keep LE arc length odd (matches A_uLE parity)
        if (k - A_uLE) % 2 != 0:
            k += 1                                            # forward = larger arc, more LE resolution
            if k >= i_LE_topo + le_offset_max:
                k -= 2
        A_lLE = k
        adapted.append('bot')

    # Top adapt: walk back from i_LE_topo - le_offset_min, find local max
    # (rare for typical turbomachinery blades; handled symmetrically).
    if margin_top < -tol:
        k = i_LE_topo - le_offset_min
        while k > i_LE_topo - le_offset_max:
            if o_t[(k - 1) % N] > o_t[k]:
                break
            k -= 1
        else:
            print(f"  WARNING: no local SS theta max within {le_offset_max} "
                  f"nodes; clamping at le_offset_max.")
        if (A_lLE - k) % 2 != 0:
            k -= 1                                            # backward = larger arc
            if k <= i_LE_topo - le_offset_max:
                k += 2
        A_uLE = k
        adapted.append('top')

    if adapted:
        print(f"  per-side LE adaptation triggered: {adapted} "
              f"(default margins: top={margin_top:+.4e}, bot={margin_bot:+.4e}, "
              f"tol={tol:.4e})")

    # Final validation
    arc = o_t[A_uLE:A_lLE + 1]
    if arc.max() > o_t[A_uLE] + tol or arc.min() < o_t[A_lLE] - tol:
        raise RuntimeError(
            f"LE arc still violates anchor band after adaptation: "
            f"arc theta range [{arc.min():.5e}, {arc.max():.5e}], "
            f"anchor band [{o_t[A_lLE]:.5e}, {o_t[A_uLE]:.5e}].")

    # Derive TE anchors so TE arc length matches LE arc length.
    le_arc_len = A_lLE - A_uLE + 1
    te_offset = (le_arc_len - 1) // 2

    # Walk te_offset outward to skip TE-tip ds-ratio kinks (see TE_DS_RATIO_TOL).
    te_offset_smooth = find_smooth_te_offset(
        o_m, o_t, te_offset, TE_DS_RATIO_TOL, TE_OFFSET_WALK_MAX, TE_KINK_BUFFER)
    if te_offset_smooth > te_offset:
        delta = te_offset_smooth - te_offset
        # Grow LE arc by 2*delta. Distribute on the side that did NOT adapt
        # during LE per-side adapt (preserves the local-extremum lock).
        if 'bot' in adapted and 'top' not in adapted:
            A_uLE -= 2 * delta                                # grow on top side
        elif 'top' in adapted and 'bot' not in adapted:
            A_lLE += 2 * delta                                # grow on bot side
        else:
            A_uLE -= delta; A_lLE += delta                    # symmetric
        print(f"  TE anchor smoothness adapt: te_offset {te_offset} -> "
              f"{te_offset_smooth}; LE arc grown by {2*delta} nodes "
              f"(A_uLE={A_uLE}, A_lLE={A_lLE}).")
        le_arc_len = A_lLE - A_uLE + 1
        te_offset = te_offset_smooth

    A_uTE = te_offset
    A_lTE = N - te_offset
    if not (0 < A_uTE < A_uLE < i_LE_topo < A_lLE < A_lTE < N):
        raise RuntimeError(
            f"anchor ordering violated: A_uTE={A_uTE}, A_uLE={A_uLE}, "
            f"i_LE_topo={i_LE_topo}, A_lLE={A_lLE}, A_lTE={A_lTE}, N={N}.")

    te_arc = np.concatenate([o_t[A_lTE:N], o_t[0:A_uTE + 1]])
    if te_arc.max() > o_t[A_uTE] + tol:
        print(f"  WARNING: TE arc theta max {te_arc.max():.5e} > A_uTE theta "
              f"{o_t[A_uTE]:.5e} (margin {o_t[A_uTE] - te_arc.max():+.3e}).")
    if te_arc.min() < o_t[A_lTE] - tol:
        print(f"  WARNING: TE arc theta min {te_arc.min():.5e} < A_lTE theta "
              f"{o_t[A_lTE]:.5e} (margin {te_arc.min() - o_t[A_lTE]:+.3e}).")

    return (i_LE_geom, i_LE_topo, A_uTE, A_uLE, A_lLE, A_lTE, N, le_arc_len)


def arclength_periodic_spline(o_m, o_t, N):
    """Periodic cubic spline on (m', theta) parametrised by arc length.
    Returns (cs_m, cs_t, s_node) with s_node[N] = total perimeter."""
    closed_m = np.r_[o_m, o_m[0]]
    closed_t = np.r_[o_t, o_t[0]]
    ds = np.hypot(np.diff(closed_m), np.diff(closed_t))
    s_node = np.concatenate(([0.0], np.cumsum(ds)))
    cs_m = CubicSpline(s_node, closed_m, bc_type='periodic')
    cs_t = CubicSpline(s_node, closed_t, bc_type='periodic')
    return cs_m, cs_t, s_node


def unit_normal_signed(dm, dt, sign):
    """Unit normal to (dm, dt). sign=+1 picks the upward (theta>0) normal,
    sign=-1 picks the downward one."""
    n1 = np.array([-dt, dm])
    n2 = np.array([dt, -dm])
    if sign > 0:
        n = n1 if n1[1] > 0 else n2
    else:
        n = n1 if n1[1] < 0 else n2
    return n / np.linalg.norm(n)


def is_node_smooth(o_m, o_t, idx, ratio_tol):
    """True when the offset-contour edge-length ratio at node idx is < ratio_tol.
    Used to pick A_uTE / A_lTE away from T-Blade3 TE-tip clustering kinks."""
    N = len(o_m)
    iL = (idx - 1) % N; iR = (idx + 1) % N
    dsL = np.hypot(o_m[idx] - o_m[iL], o_t[idx] - o_t[iL])
    dsR = np.hypot(o_m[iR]  - o_m[idx], o_t[iR]  - o_t[idx])
    if dsL == 0.0 or dsR == 0.0:
        return False
    return max(dsL, dsR) / min(dsL, dsR) < ratio_tol


def find_smooth_te_offset(o_m, o_t, te_offset_initial, ratio_tol, walk_max,
                          buffer):
    """Walk te_offset outward (parity-preserving, +2 per step) until A_uTE,
    A_lTE AND `buffer` nodes interior to each (toward the closure) are all
    smooth. Returns the new te_offset, or te_offset_initial if no valid
    position is found within walk_max steps."""
    N = len(o_m)
    for te in range(te_offset_initial, te_offset_initial + walk_max + 1, 2):
        all_smooth = True
        for k in range(buffer + 1):
            if not is_node_smooth(o_m, o_t, te - k, ratio_tol):
                all_smooth = False; break
            if not is_node_smooth(o_m, o_t, N - te + k, ratio_tol):
                all_smooth = False; break
        if all_smooth:
            return te
    print(f"  WARNING: no smooth TE-anchor pair (buffer={buffer}) within "
          f"walk_max={walk_max} of te_offset_initial={te_offset_initial}; "
          f"falling back.")
    return te_offset_initial


def inflate_arc(cs_m, cs_t, s_node, i_start, i_end, dist, sign):
    """Pointwise offset of contour nodes [i_start..i_end] by dist along the
    local signed normal. Returns (n_pts, 2)."""
    n_pts = i_end - i_start + 1
    out = np.zeros((n_pts, 2))
    for k, i_local in enumerate(range(i_start, i_end + 1)):
        s_k = s_node[i_local]
        dm_k = float(cs_m(s_k, 1))
        dt_k = float(cs_t(s_k, 1))
        n_k = unit_normal_signed(dm_k, dt_k, sign)
        out[k, 0] = float(cs_m(s_k)) + dist * n_k[0]
        out[k, 1] = float(cs_t(s_k)) + dist * n_k[1]
    return out


def one_sided_tanh(N, L, ds_target, which):
    """N points on [0, L] with a tanh distribution clustered at i=N-1
    (last segment is the smallest, first segment is the largest).
    `which` specifies which segment is constrained to ds_target:
        'last'  -> last segment  (i=N-2 -> i=N-1, the SMALLEST). Used by
                   the L-block horizontal clustering: ds_target is the
                   spacing approaching the LE/TE corner.
        'first' -> first segment (i=0   -> i=1,   the LARGEST). Used by
                   the O-grid radial cut: ds_target is the outer-side
                   first step (eta=0 -> eta=1).
    To cluster at i=0 instead, call with reversed parameters and post-flip
    the result (1 - x[::-1])."""
    u = np.arange(N) / (N - 1)

    def positions(tau):
        return L * np.tanh(tau * u) / np.tanh(tau)

    if which == 'last':
        if ds_target >= L / (N - 1) - 1e-12:
            return np.linspace(0.0, L, N)
        def constraint(tau):
            x = positions(tau)
            return (x[-1] - x[-2]) - ds_target
    elif which == 'first':
        # Constraining the largest segment: meaningful only when
        # ds_target > L/(N-1) (otherwise the clustering is the wrong sign).
        if ds_target <= L / (N - 1) + 1e-12:
            return np.linspace(0.0, L, N)
        def constraint(tau):
            x = positions(tau)
            return (x[1] - x[0]) - ds_target
    else:
        raise ValueError("which must be 'first' or 'last'")
    tau = brentq(constraint, 1e-6, 50.0)
    return positions(tau)


# ============================================================================
# I/O
# ============================================================================

CASE_NEWTON_TEMPLATE = """\
# ============================================================================
# gridexec input -- {label} ({nx} x {ny})
# ============================================================================

mode             = 2D
niter            = {niter}
nblocks          = 1
npts_x           = {nx}
npts_y           = {ny}

solver           = newton

lambda_max       = {lam}
ramp_iters       = {ramp}

tol_step         = {tol_step}
tol_resid        = {tol_resid}
"""

CASE_PSN_LIM_TEMPLATE = """\
# ============================================================================
# gridexec input -- {label} ({nx} x {ny})
# ============================================================================

mode             = 2D
niter            = {niter}
nblocks          = 1
npts_x           = {nx}
npts_y           = {ny}

solver           = psn_lim

lambda_max       = {lam}
ramp_iters       = {ramp}
inner_niter      = {inner}

tol_step         = {tol_step}
tol_resid        = {tol_resid}
"""

FMT_4COL = "%24.16e %24.16e %24.16e %24.16e"


def write_block(out_dir, label, bot, top, left, right, case_text):
    """Write x_boundaries.dat, y_boundaries.dat, case.in for one block.
    bot, top : (nx, 2). left, right : (ny, 2). Corner consistency NOT
    re-checked here — caller is expected to have validated."""
    nx = len(bot)
    ny = len(left)
    assert len(top) == nx, f"{label}: top has {len(top)} pts, expected {nx}"
    assert len(right) == ny, f"{label}: right has {len(right)} pts, expected {ny}"
    os.makedirs(out_dir, exist_ok=True)
    np.savetxt(os.path.join(out_dir, "x_boundaries.dat"),
               np.column_stack([bot, top]), fmt=FMT_4COL)
    np.savetxt(os.path.join(out_dir, "y_boundaries.dat"),
               np.column_stack([left, right]), fmt=FMT_4COL)
    with open(os.path.join(out_dir, "case.in"), "w") as f:
        f.write(case_text)


def check_corners(label, bot, top, left, right, tol=1e-12):
    pairs = [
        ("(eta=0,xi=0)   bot[0]   == left[0]   ",  bot[0],   left[0]),
        ("(eta=0,xi=N-1) bot[-1]  == right[0]  ",  bot[-1],  right[0]),
        ("(eta=1,xi=0)   top[0]   == left[-1]  ",  top[0],   left[-1]),
        ("(eta=1,xi=N-1) top[-1]  == right[-1] ",  top[-1],  right[-1]),
    ]
    print(f"\n[{label}] corner consistency:")
    for desc, p, q in pairs:
        d = float(np.linalg.norm(p - q))
        flag = "OK" if d < tol else "MISMATCH"
        print(f"  {desc}: |p-q|={d:.3e}  {flag}")


# ============================================================================
# Per-block builders
# ============================================================================

def build_ogrid(blade_closed, ogrid_closed):
    """O-grid block (single cut). xi=0 and xi=1 walls are identical, going
    from ogrid[0] (eta=0, outer, TE) to blade[0] (eta=1, inner, TE)."""
    assert blade_closed.shape == ogrid_closed.shape
    n_xi = blade_closed.shape[0]                              # closed count
    n_eta = N_OGRID_RADIAL

    P_outer = ogrid_closed[0].copy()
    P_inner = blade_closed[0].copy()
    L_cut = float(np.linalg.norm(P_inner - P_outer))

    # ogrid clustering: cluster at eta=N-1 (blade), specify the FIRST
    # (outer-side) segment as DS_OGRID_FIRST. This was the working OGV55
    # convention — see project_ogrid_3.md.
    pos = one_sided_tanh(n_eta, L_cut, DS_OGRID_FIRST, which='first')
    u_cut = pos / L_cut
    cut = P_outer + np.outer(u_cut, P_inner - P_outer)        # (n_eta, 2)

    # eta=0 row = outer (ogrid), eta=1 row = inner (blade); single cut -> both
    # xi walls are the same.
    bot = ogrid_closed.copy()
    top = blade_closed.copy()
    left = cut.copy()
    right = cut.copy()
    return bot, top, left, right, n_xi, n_eta, L_cut


def build_passage_blocks(ogrid_closed, A_uTE, A_uLE, A_lLE, A_lTE, N, n_y_passage):
    """Build the inflated lids and shared L-block clustering distributions
    used by both top and bot passages. Returns dict with all derived arrays
    and scalars needed by the per-passage builders."""
    o_m = ogrid_closed[:N, 0]
    o_t = ogrid_closed[:N, 1]
    cs_m, cs_t, s_node = arclength_periodic_spline(o_m, o_t, N)

    # Anchor coordinates on the offset contour
    s_uLE = s_node[A_uLE]; s_uTE = s_node[A_uTE]
    s_lLE = s_node[A_lLE]; s_lTE = s_node[A_lTE]
    m_uLE, th_uLE = float(cs_m(s_uLE)), float(cs_t(s_uLE))
    m_uTE, th_uTE = float(cs_m(s_uTE)), float(cs_t(s_uTE))
    m_lLE, th_lLE = float(cs_m(s_lLE)), float(cs_t(s_lLE))
    m_lTE, th_lTE = float(cs_m(s_lTE)), float(cs_t(s_lTE))

    # Inflated lids (pointwise offset along local normal)
    inflated_top = inflate_arc(cs_m, cs_t, s_node, A_uTE, A_uLE,
                               INFLATE_DIST_TOP, sign=+1)
    inflated_bot = inflate_arc(cs_m, cs_t, s_node, A_lLE, A_lTE,
                               INFLATE_DIST_BOT, sign=-1)

    # Lid endpoints
    m_lid_LE_top, th_lid_LE_top = float(inflated_top[-1, 0]), float(inflated_top[-1, 1])
    m_lid_TE_top, th_lid_TE_top = float(inflated_top[ 0, 0]), float(inflated_top[ 0, 1])
    m_lid_LE_bot, th_lid_LE_bot = float(inflated_bot[ 0, 0]), float(inflated_bot[ 0, 1])
    m_lid_TE_bot, th_lid_TE_bot = float(inflated_bot[-1, 0]), float(inflated_bot[-1, 1])

    # Reference SS / PS arcs for the inflated-lid coords used downstream
    upperSS = ogrid_closed[A_uTE:A_uLE + 1].copy()            # A_uTE -> A_uLE
    lowerPS = ogrid_closed[A_lLE:A_lTE + 1].copy()            # A_lLE -> A_lTE

    # Shared L-block lengths (averaged top + bot lid distance to the wall, same
    # rule as the original build scripts so the LEFT/RIGHT L-blocks of top
    # and bot passage share their normalised i-distribution).
    L_LE_top = m_lid_LE_top - M_LEFT_END
    L_LE_bot = m_lid_LE_bot - M_LEFT_END
    L_TE_top = M_RIGHT_END - m_lid_TE_top
    L_TE_bot = M_RIGHT_END - m_lid_TE_bot
    L_LE = 0.5 * (L_LE_top + L_LE_bot)
    L_TE = 0.5 * (L_TE_top + L_TE_bot)

    pos_ref_LE = one_sided_tanh(N_X_LR, L_LE, SHARED_LE_DS, which='last')
    pos_ref_TE = one_sided_tanh(N_X_LR, L_TE, SHARED_TE_DS, which='last')
    u_LE = pos_ref_LE / L_LE                                  # cluster at i=N-1
    u_TE = pos_ref_TE / L_TE                                  # cluster at i=N-1

    return dict(
        m_uLE=m_uLE, th_uLE=th_uLE, m_uTE=m_uTE, th_uTE=th_uTE,
        m_lLE=m_lLE, th_lLE=th_lLE, m_lTE=m_lTE, th_lTE=th_lTE,
        upperSS=upperSS, lowerPS=lowerPS,
        inflated_top=inflated_top, inflated_bot=inflated_bot,
        m_lid_LE_top=m_lid_LE_top, th_lid_LE_top=th_lid_LE_top,
        m_lid_TE_top=m_lid_TE_top, th_lid_TE_top=th_lid_TE_top,
        m_lid_LE_bot=m_lid_LE_bot, th_lid_LE_bot=th_lid_LE_bot,
        m_lid_TE_bot=m_lid_TE_bot, th_lid_TE_bot=th_lid_TE_bot,
        u_LE=u_LE, u_TE=u_TE, L_LE=L_LE, L_TE=L_TE,
        n_y_passage=n_y_passage,
    )


def build_top_passage(d):
    """Top passage: 3 blocks (LEFT, CENTER, RIGHT). eta=0 = blade-anchor row,
    eta=1 = inflated upperSS lid above."""
    n_y_passage = d['n_y_passage']
    u_LE = d['u_LE']; u_TE = d['u_TE']
    th_LE = d['th_uLE']; th_TE = d['th_uTE']
    th_lid_LE = d['th_lid_LE_top']; th_lid_TE = d['th_lid_TE_top']
    m_LE = d['m_uLE']; m_TE = d['m_uTE']
    m_lid_LE = d['m_lid_LE_top']; m_lid_TE = d['m_lid_TE_top']

    # ----- LEFT -----
    L_top_L = m_lid_LE - M_LEFT_END
    L_bot_L = m_LE     - M_LEFT_END
    left_top  = np.column_stack([M_LEFT_END + u_LE * L_top_L,
                                 np.full(N_X_LR, th_lid_LE)])
    left_bot  = np.column_stack([M_LEFT_END + u_LE * L_bot_L,
                                 np.full(N_X_LR, th_LE)])
    left_wall = np.column_stack([np.full(n_y_passage, M_LEFT_END),
                                 np.linspace(th_LE, th_lid_LE, n_y_passage)])
    red_curve = np.column_stack([np.linspace(m_LE, m_lid_LE, n_y_passage),
                                 np.linspace(th_LE, th_lid_LE, n_y_passage)])

    # ----- CENTER -----
    # upperSS goes A_uTE -> A_uLE; reverse so xi=0 is at the LE side, matching
    # LEFT.right and RIGHT.left.
    center_bot  = d['upperSS'][::-1].copy()
    center_top  = d['inflated_top'][::-1].copy()
    center_left = red_curve.copy()
    blue_curve  = np.column_stack([np.linspace(m_TE, m_lid_TE, n_y_passage),
                                   np.linspace(th_TE, th_lid_TE, n_y_passage)])

    # ----- RIGHT -----
    u_TE_rev = 1.0 - u_TE[::-1]                               # cluster at i=0 (TE corner)
    L_top_R = M_RIGHT_END - m_lid_TE
    L_bot_R = M_RIGHT_END - m_TE
    right_top  = np.column_stack([m_lid_TE + u_TE_rev * L_top_R,
                                  np.full(N_X_LR, th_lid_TE)])
    right_bot  = np.column_stack([m_TE     + u_TE_rev * L_bot_R,
                                  np.full(N_X_LR, th_TE)])
    right_left = blue_curve.copy()
    right_wall = np.column_stack([np.full(n_y_passage, M_RIGHT_END),
                                  np.linspace(th_TE, th_lid_TE, n_y_passage)])

    return {
        "left":   (left_bot,   left_top,   left_wall,   red_curve),
        "center": (center_bot, center_top, center_left, blue_curve),
        "right":  (right_bot,  right_top,  right_left,  right_wall),
    }


def build_bot_passage(d):
    """Bot passage: 3 blocks. eta=0 = inflated-down lid (lower theta),
    eta=1 = lowerPS-anchor row (higher theta)."""
    n_y_passage = d['n_y_passage']
    u_LE = d['u_LE']; u_TE = d['u_TE']
    th_lLE = d['th_lLE']; th_lTE = d['th_lTE']
    th_lid_LE = d['th_lid_LE_bot']; th_lid_TE = d['th_lid_TE_bot']
    m_lLE = d['m_lLE']; m_lTE = d['m_lTE']
    m_lid_LE = d['m_lid_LE_bot']; m_lid_TE = d['m_lid_TE_bot']

    # ----- LEFT -----
    L_bot_L = m_lid_LE - M_LEFT_END
    L_top_L = m_lLE    - M_LEFT_END
    left_bot  = np.column_stack([M_LEFT_END + u_LE * L_bot_L,
                                 np.full(N_X_LR, th_lid_LE)])
    left_top  = np.column_stack([M_LEFT_END + u_LE * L_top_L,
                                 np.full(N_X_LR, th_lLE)])
    left_wall = np.column_stack([np.full(n_y_passage, M_LEFT_END),
                                 np.linspace(th_lid_LE, th_lLE, n_y_passage)])
    red_curve = np.column_stack([np.linspace(m_lid_LE, m_lLE, n_y_passage),
                                 np.linspace(th_lid_LE, th_lLE, n_y_passage)])

    # ----- CENTER -----
    # inflated_bot already runs A_lLE -> A_lTE; lowerPS too.
    center_bot  = d['inflated_bot'].copy()
    center_top  = d['lowerPS'].copy()
    center_left = red_curve.copy()
    blue_curve  = np.column_stack([np.linspace(m_lid_TE, m_lTE, n_y_passage),
                                   np.linspace(th_lid_TE, th_lTE, n_y_passage)])

    # ----- RIGHT -----
    u_TE_rev = 1.0 - u_TE[::-1]
    L_bot_R = M_RIGHT_END - m_lid_TE
    L_top_R = M_RIGHT_END - m_lTE
    right_bot  = np.column_stack([m_lid_TE + u_TE_rev * L_bot_R,
                                  np.full(N_X_LR, th_lid_TE)])
    right_top  = np.column_stack([m_lTE    + u_TE_rev * L_top_R,
                                  np.full(N_X_LR, th_lTE)])
    right_left = blue_curve.copy()
    right_wall = np.column_stack([np.full(n_y_passage, M_RIGHT_END),
                                  np.linspace(th_lid_TE, th_lTE, n_y_passage)])

    return {
        "left":   (left_bot,   left_top,   left_wall,   red_curve),
        "center": (center_bot, center_top, center_left, blue_curve),
        "right":  (right_bot,  right_top,  right_left,  right_wall),
    }


def build_middle_passage(ogrid_closed, A_uTE, A_uLE, A_lLE, A_lTE, N,
                         n_y_passage, top_blocks, bot_blocks):
    """Middle blocks closing the gap between top and bot at the LE and TE.
    Reuses the bot.left.top + top.left.bot horizontals (and similarly for
    right) so eta-lines are not slanted."""
    # ----- middle.left -----
    # eta=0 = bot.left.top  (bot passage's TOP horizontal of LEFT block)
    # eta=1 = top.left.bot  (top passage's BOT horizontal of LEFT block)
    ml_bot = bot_blocks["left"][1].copy()                     # bot.left.top
    ml_top = top_blocks["left"][0].copy()                     # top.left.bot
    # right wall (xi=1) = LE arc, traversed A_lLE -> A_uLE (so eta increases
    # with theta: lower-theta endpoint at eta=0)
    ml_right = ogrid_closed[A_uLE:A_lLE + 1][::-1].copy()
    if ml_right.shape[0] != n_y_passage:
        raise RuntimeError(
            f"middle.left right wall has {ml_right.shape[0]} pts, "
            f"expected {n_y_passage}. LE arc length / n_y_passage are "
            f"derived from find_anchors — internal inconsistency.")
    # left wall (xi=0) = vertical at M_LEFT_END
    th_lo = float(ml_bot[0, 1])
    th_hi = float(ml_top[0, 1])
    ml_left = np.column_stack([np.full(n_y_passage, M_LEFT_END),
                               np.linspace(th_lo, th_hi, n_y_passage)])

    # ----- middle.right -----
    mr_bot = bot_blocks["right"][1].copy()                    # bot.right.top
    mr_top = top_blocks["right"][0].copy()                    # top.right.bot
    # left wall (xi=0) = TE arc, traversed A_lTE -> A_uTE through the wrap
    mr_left = np.r_[ogrid_closed[A_lTE:N], ogrid_closed[0:A_uTE + 1]]
    if mr_left.shape[0] != n_y_passage:
        raise RuntimeError(
            f"middle.right left wall has {mr_left.shape[0]} pts, "
            f"expected {n_y_passage}. TE arc length and n_y_passage are "
            f"derived from find_anchors — internal inconsistency.")
    th_lo = float(mr_bot[-1, 1])
    th_hi = float(mr_top[-1, 1])
    mr_right = np.column_stack([np.full(n_y_passage, M_RIGHT_END),
                                np.linspace(th_lo, th_hi, n_y_passage)])

    return {
        "left":  (ml_bot, ml_top, ml_left, ml_right),
        "right": (mr_bot, mr_top, mr_left, mr_right),
    }


# ============================================================================
# Driver
# ============================================================================

def _format_tol(x):
    """Format a tolerance like Python would but without zero-padded exponent
    (matches the hand-written originals: 1.0e-8 instead of 1.0e-08)."""
    s = f"{x:.1e}"                                            # e.g. '1.0e-08'
    mant, exp = s.split('e')
    return f"{mant}e{int(exp)}"                               # '1.0e-8'


def newton_case_text(label, nx, ny, lam):
    return CASE_NEWTON_TEMPLATE.format(
        label=label, nx=nx, ny=ny, niter=NITER_NEWTON,
        lam=lam, ramp=RAMP_ITERS_NEWTON,
        tol_step=_format_tol(TOL_STEP), tol_resid=_format_tol(TOL_RESID),
    )


def psn_lim_case_text(label, nx, ny, lam):
    return CASE_PSN_LIM_TEMPLATE.format(
        label=label, nx=nx, ny=ny, niter=NITER_OGRID,
        lam=lam, ramp=RAMP_ITERS_OGRID, inner=INNER_NITER_OGRID,
        tol_step=_format_tol(TOL_STEP), tol_resid=_format_tol(TOL_RESID),
    )


def plot_topology_preview(out_path, blade_closed, ogrid_closed, derived,
                          ogrid_block, top_blocks, bot_blocks, mid_blocks,
                          dpi=140):
    """Pre-solve sanity check. Draws blade + offset contour as context, the
    inflated lids, and the 4 walls of every block in distinct colors per
    passage. Runs after all 9 blocks have been built; before any solver."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(20, 9))

    # context
    ax.plot(blade_closed[:, 0], blade_closed[:, 1],
            color='gray', lw=1.0, alpha=0.6, label='blade')
    ax.plot(ogrid_closed[:, 0], ogrid_closed[:, 1],
            color='gray', lw=0.5, ls='--', alpha=0.5,
            label='ogrid contour (eta=0 of O-grid block)')
    ax.plot(derived['inflated_top'][:, 0], derived['inflated_top'][:, 1],
            color='C1', lw=0.8, ls=':', alpha=0.7, label='inflated lid (top)')
    ax.plot(derived['inflated_bot'][:, 0], derived['inflated_bot'][:, 1],
            color='C0', lw=0.8, ls=':', alpha=0.7, label='inflated lid (bot)')

    def draw_block(walls, color, label=None):
        bot, top, left, right = walls
        for i, w in enumerate((bot, top, left, right)):
            ax.plot(w[:, 0], w[:, 1], color=color, lw=0.6,
                    label=label if i == 0 else None)

    # ogrid radial cut (the only non-context wall of the O-grid block)
    o_left = ogrid_block[2]
    ax.plot(o_left[:, 0], o_left[:, 1], color='C2', lw=0.8,
            label='ogrid_3 radial cut')

    for i, name in enumerate(('left', 'center', 'right')):
        draw_block(top_blocks[name], 'C1',
                   label='top_passage_3' if i == 0 else None)
        draw_block(bot_blocks[name], 'C0',
                   label='bot_passage_3' if i == 0 else None)
    for i, name in enumerate(('left', 'right')):
        draw_block(mid_blocks[name], 'C3',
                   label='middle_passage_3' if i == 0 else None)

    ax.set_aspect('equal')
    ax.grid(alpha=0.3)
    ax.set_xlabel("m'")
    ax.set_ylabel(r"$\theta$")
    ax.set_title("9-block topology preview (boundaries only — pre-solve)")
    ax.legend(loc='best', fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close(fig)


def write_passage(out_root, passage_name, blocks, lambda_map):
    section = os.path.basename(os.path.normpath(out_root))
    for block_name in ("left", "center", "right"):
        if block_name not in blocks:
            continue
        bot, top, left, right = blocks[block_name]
        nx, ny = len(bot), len(left)
        check_corners(f"{passage_name}.{block_name}", bot, top, left, right)
        label = f"{section}/{passage_name}/{block_name}"
        case_text = newton_case_text(
            label=label, nx=nx, ny=ny, lam=lambda_map[block_name])
        write_block(os.path.join(out_root, passage_name, block_name),
                    label=label, bot=bot, top=top, left=left, right=right,
                    case_text=case_text)
        print(f"wrote {passage_name}/{block_name}: "
              f"{nx}x{ny}, lambda_max={lambda_map[block_name]}")


def main():
    global INFLATE_DIST_TOP, INFLATE_DIST_BOT
    parser = argparse.ArgumentParser(
        description="Build the 9-block blade-section grid topology.")
    parser.add_argument("blade_path",
                        help="Path to a T-Blade3 blade.X.X.SECTION file")
    parser.add_argument("output_dir",
                        help="Root directory where block subdirs are written")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip the topology_preview.png sanity-check plot")
    parser.add_argument("--mode", choices=["adaptive", "unified"],
                        default="adaptive",
                        help="adaptive (default): per-section LE+TE adaptation. "
                             "unified: freeze anchors via --unified-anchors so "
                             "all sections produce stackable identical-dim blocks.")
    parser.add_argument("--unified-anchors", metavar="A_uTE,A_uLE,A_lLE,A_lTE",
                        help="Required with --mode unified. Comma-separated ints. "
                             "For OGV55: 20,88,128,220.")
    parser.add_argument("--inflate-dist", type=float, default=None,
                        help="Override INFLATE_DIST_TOP/BOT. Recommended for "
                             "unified mode where 0.10 default folds bot.center "
                             "(use 0.06 for OGV55).")
    args = parser.parse_args()

    if args.inflate_dist is not None:
        INFLATE_DIST_TOP = args.inflate_dist
        INFLATE_DIST_BOT = args.inflate_dist

    # ----- T-Blade3 input -----
    m, t, chord, pitch = read_tblade3_blade(args.blade_path)
    delta = OFFSET_FRAC * chord
    print(f"T-Blade3 file: {args.blade_path}")
    print(f"  chord = {chord:.6f}, pitch = {pitch:.6f}, "
          f"unique blade pts = {len(m)}, offset = {delta:.6e}")
    print(f"  mode = {args.mode}, "
          f"INFLATE_DIST_TOP/BOT = {INFLATE_DIST_TOP}/{INFLATE_DIST_BOT}")

    # Reverse blade if it came in CW (we need CCW for the offset/anchor logic).
    blade_closed = np.column_stack([np.r_[m, m[0]], np.r_[t, t[0]]])
    if signed_area(blade_closed) <= 0:
        print("  blade contour was CW; reversing to CCW")
        m = m[::-1]; t = t[::-1]
        blade_closed = np.column_stack([np.r_[m, m[0]], np.r_[t, t[0]]])

    # ----- O-grid contour -----
    ogrid_closed = build_ogrid_contour(m, t, delta)
    assert signed_area(ogrid_closed) > 0

    # ----- Anchors (adaptive or unified) -----
    if args.mode == "unified":
        if args.unified_anchors is None:
            parser.error("--mode unified requires --unified-anchors A_uTE,A_uLE,A_lLE,A_lTE")
        try:
            vals = [int(x) for x in args.unified_anchors.split(",")]
        except ValueError:
            parser.error("--unified-anchors must be 4 comma-separated integers")
        if len(vals) != 4:
            parser.error("--unified-anchors must have 4 values")
        A_uTE, A_uLE, A_lLE, A_lTE = vals
        N = ogrid_closed.shape[0] - 1
        i_LE_topo = N // 2
        i_LE_geom = int(ogrid_closed[:N, 0].argmin())
        # Validate: ordering, theta-band, TE smoothness with buffer
        if not (0 < A_uTE < A_uLE < i_LE_topo < A_lLE < A_lTE < N):
            raise RuntimeError(
                f"unified anchors fail ordering for N={N}: "
                f"A_uTE={A_uTE}, A_uLE={A_uLE}, i_LE_topo={i_LE_topo}, "
                f"A_lLE={A_lLE}, A_lTE={A_lTE}")
        if (A_lLE - A_uLE + 1) != (2 * A_uTE + 1):
            raise RuntimeError(
                f"unified LE arc length ({A_lLE-A_uLE+1}) must equal "
                f"2*A_uTE+1 ({2*A_uTE+1}) for n_y_passage to match between "
                f"middle.left (LE arc) and middle.right (TE arc).")
        n_y_passage = A_lLE - A_uLE + 1
        o_m = ogrid_closed[:N, 0]; o_t = ogrid_closed[:N, 1]
        tol = LE_INTERSECTION_TOL_FRAC * chord
        arc = o_t[A_uLE:A_lLE + 1]
        margin_top = o_t[A_uLE] - arc.max()
        margin_bot = arc.min() - o_t[A_lLE]
        if margin_top < -tol or margin_bot < -tol:
            raise RuntimeError(
                f"unified anchors fail theta-band for this section: "
                f"margin_top={margin_top:+.3e}, margin_bot={margin_bot:+.3e}, "
                f"tol={tol:.3e}. The anchors don't bound this section's LE arc.")
        te_smooth = all(
            is_node_smooth(o_m, o_t, A_uTE - k, TE_DS_RATIO_TOL)
            and is_node_smooth(o_m, o_t, A_lTE + k, TE_DS_RATIO_TOL)
            for k in range(TE_KINK_BUFFER + 1))
        if not te_smooth:
            print(f"  WARNING: TE smoothness with buffer={TE_KINK_BUFFER} "
                  f"fails at unified A_uTE/A_lTE for this section; centers "
                  f"may fold. Consider per-section adaptive mode.")
    else:
        (i_LE_geom, i_LE_topo, A_uTE, A_uLE, A_lLE, A_lTE, N,
         n_y_passage) = find_anchors(
            ogrid_closed, chord, LE_OFFSET_DEFAULT, LE_OFFSET_MIN,
            LE_OFFSET_MAX, LE_INTERSECTION_TOL_FRAC, LE_GEOM_TOPO_TOL)
    print(f"  N={N}, i_LE_geom = {i_LE_geom}, i_LE_topo = {i_LE_topo}, "
          f"A_uTE = {A_uTE}, A_uLE = {A_uLE}, "
          f"A_lLE = {A_lLE}, A_lTE = {A_lTE}")
    print(f"  LE arc = TE arc = n_y_passage = {n_y_passage} pts")

    # ----- Passage-shared geometry (lids + L-block distributions) -----
    derived = build_passage_blocks(ogrid_closed, A_uTE, A_uLE, A_lLE, A_lTE,
                                   N, n_y_passage)
    print(f"  L_LE = {derived['L_LE']:.6f}  (avg top+bot lid -> M_LEFT_END), "
          f"SHARED_LE_DS = {SHARED_LE_DS:.6e}")
    print(f"  L_TE = {derived['L_TE']:.6f}  (avg top+bot lid -> M_RIGHT_END), "
          f"SHARED_TE_DS = {SHARED_TE_DS:.6e}")

    # ----- Build all 9 blocks -----
    top_blocks = build_top_passage(derived)
    bot_blocks = build_bot_passage(derived)
    mid_blocks = build_middle_passage(
        ogrid_closed, A_uTE, A_uLE, A_lLE, A_lTE, N, n_y_passage,
        top_blocks, bot_blocks)

    # ----- Write -----
    out_root = os.path.expanduser(args.output_dir)
    os.makedirs(out_root, exist_ok=True)

    # ogrid_3
    bot, top, left, right, n_xi, n_eta, L_cut = build_ogrid(
        blade_closed, ogrid_closed)
    print(f"\nogrid_3: radial cut L = {L_cut:.6f}, "
          f"first ds = {DS_OGRID_FIRST}")
    check_corners("ogrid_3", bot, top, left, right)
    section = os.path.basename(os.path.normpath(out_root))
    label = f"{section}/ogrid_3"
    write_block(os.path.join(out_root, "ogrid_3"),
                label=label, bot=bot, top=top, left=left, right=right,
                case_text=psn_lim_case_text(
                    label=label, nx=n_xi, ny=n_eta, lam=LAMBDA_OGRID))
    print(f"wrote ogrid_3: {n_xi}x{n_eta} (psn_lim, lambda_max={LAMBDA_OGRID})")

    # passages
    print()
    write_passage(out_root, "top_passage_3", top_blocks, LAMBDA_TOP)
    write_passage(out_root, "bot_passage_3", bot_blocks, LAMBDA_BOT)
    write_passage(out_root, "middle_passage_3", mid_blocks, LAMBDA_MID)

    # preview plot
    if not args.no_plot:
        ogrid_block = (bot, top, left, right)
        preview_path = os.path.join(out_root, "topology_preview.png")
        plot_topology_preview(preview_path, blade_closed, ogrid_closed,
                              derived, ogrid_block,
                              top_blocks, bot_blocks, mid_blocks)
        print(f"\nwrote topology preview: {preview_path}")

    print(f"\nDone. Output root: {out_root}")


if __name__ == "__main__":
    main()
