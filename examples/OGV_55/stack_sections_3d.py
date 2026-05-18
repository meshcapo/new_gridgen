"""Stack 21 (m', theta) blade-section grids into a 3D 9-block CFD grid by
applying T-Blade3's (m', theta) -> (x, y, z) transformation per section.

Per section ia (1..nsl):
    m'_3D = m'_blade + dmp(ia) + delmp(ia)        # delmp = sweep
    theta_3D = theta_blade + delta_theta(ia)      # lean
    x = x_streamline_ia(m'_3D)
    r = r_streamline_ia(m'_3D)
    y = r * sin(theta_3D)
    z = r * cos(theta_3D)

dmp(ia) = msle(ia) - m'_blade(LE, ia), where msle is the m'-coord on
streamline ia at the intersection with the LE line (defined by the LE/TE
input curve in 3dbgbinput).

Each of the 9 blocks becomes a 3D structured grid (.vts) with shape
(nx, ny, nsl). PointData fields (orthogonality, aspect_ratio, jacobian) from
each section are stacked along the spanwise axis.

Usage:
    python3 stack_sections_3d.py <sections_root> <3dbgbinput.path> <output_dir>
        [--vts-name final_grid.vts] [--le-idx 120]

Writes <output_dir>/stacked_3d.vtm referencing 9 stacked .vts files.
"""
import os
import sys
import argparse
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from vtkmodules.vtkIOXML import (vtkXMLStructuredGridReader,
                                  vtkXMLStructuredGridWriter)
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid
from vtkmodules.vtkCommonCore import vtkPoints, vtkFloatArray
from vtkmodules.util.numpy_support import (vtk_to_numpy, numpy_to_vtk)


BLOCKS = [
    ("ogrid_3",         "ogrid_3"),
    ("top.left",        "top_passage_3/left"),
    ("top.center",      "top_passage_3/center"),
    ("top.right",       "top_passage_3/right"),
    ("bot.left",        "bot_passage_3/left"),
    ("bot.center",      "bot_passage_3/center"),
    ("bot.right",       "bot_passage_3/right"),
    ("middle.left",     "middle_passage_3/left"),
    ("middle.right",    "middle_passage_3/right"),
]


def parse_3dbgbinput(path):
    """Return dict with nsl, sweep_cps [(span,delta_m)], lean_cps,
    LE/TE line endpoints, and a list of nsl streamlines (each (npts, 2))."""
    with open(path) as f:
        lines = f.readlines()

    def find(keyword, start=0):
        for i in range(start, len(lines)):
            if keyword in lines[i]:
                return i
        raise KeyError(f"keyword not found: {keyword}")

    # nsl
    i = find("Number of streamlines")
    nsl = int(lines[i + 1].split()[0])

    # LE/TE curve endpoints (assume 2-point linear definition)
    i = find("LE / TE curve")
    n_curve = int(lines[i + 2].split()[0])
    if n_curve != 2:
        raise NotImplementedError(
            f"LE/TE curve with {n_curve} points; only 2-point linear supported "
            f"(extend by interpolation if needed).")
    hub_vals = [float(x) for x in lines[i + 4].split()[:4]]
    cas_vals = [float(x) for x in lines[i + 5].split()[:4]]
    le_hub = (hub_vals[0], hub_vals[1])
    le_cas = (cas_vals[0], cas_vals[1])
    te_hub = (hub_vals[2], hub_vals[3])
    te_cas = (cas_vals[2], cas_vals[3])

    # Sweep
    i = find("Control points for sweep")
    n_sw = int(lines[i + 1].split()[0])
    sweep_cps = np.array([[float(x) for x in lines[i + 3 + k].split()[:2]]
                          for k in range(n_sw)])

    # Lean
    i = find("Control points for lean")
    n_ln = int(lines[i + 1].split()[0])
    lean_cps = np.array([[float(x) for x in lines[i + 3 + k].split()[:2]]
                         for k in range(n_ln)])

    # Streamlines: skip until we see "x_s" header, then read nsl blocks
    # delimited by "0 0" sentinels
    i = find("Streamline Data")
    while not lines[i].strip().startswith("x_s"):
        i += 1
    i += 1                                                        # past header
    streams = []
    cur = []
    while i < len(lines) and len(streams) < nsl:
        parts = lines[i].split()
        if len(parts) >= 2:
            try:
                x, r = float(parts[0]), float(parts[1])
            except ValueError:
                break
            if x == 0.0 and r == 0.0:
                if cur:
                    streams.append(np.array(cur)); cur = []
            else:
                cur.append((x, r))
        elif not parts:
            pass
        else:
            break
        i += 1
    if cur and len(streams) < nsl:
        streams.append(np.array(cur))
    if len(streams) != nsl:
        raise RuntimeError(
            f"expected {nsl} streamlines, parsed {len(streams)}.")
    return {
        "nsl": nsl, "sweep_cps": sweep_cps, "lean_cps": lean_cps,
        "le_hub": le_hub, "le_cas": le_cas,
        "te_hub": te_hub, "te_cas": te_cas, "streams": streams,
    }


def streamline_splines(streamline):
    """Return (mp_s, x_spl, r_spl). m' is computed per T-Blade3:
    dm = 2*sqrt(dx^2 + dr^2) / (r_i + r_{i-1})."""
    x_s, r_s = streamline[:, 0], streamline[:, 1]
    mp = np.zeros(len(streamline))
    for k in range(1, len(streamline)):
        mp[k] = mp[k - 1] + (2.0 * np.hypot(x_s[k] - x_s[k - 1],
                                             r_s[k] - r_s[k - 1])
                             / (r_s[k] + r_s[k - 1]))
    return mp, CubicSpline(mp, x_s), CubicSpline(mp, r_s)


def line_intersection_msle(mp_s, x_spl, r_spl, hub, cas):
    """Find m' on streamline where it crosses the line through hub..cas.
    Cross-product sign-change gives the intersection."""
    dx = cas[0] - hub[0]; dr = cas[1] - hub[1]

    def f(m):
        return (x_spl(m) - hub[0]) * dr - (r_spl(m) - hub[1]) * dx

    for k in range(len(mp_s) - 1):
        if f(mp_s[k]) * f(mp_s[k + 1]) < 0:
            return brentq(f, mp_s[k], mp_s[k + 1])
    raise RuntimeError("streamline does not cross the LE/TE line in its m' range.")


def interp_cps(span, cps):
    """Linear interpolation of (span, value) control points. Falls back to
    endpoint values outside the table."""
    spans, vals = cps[:, 0], cps[:, 1]
    return float(np.interp(span, spans, vals))


def read_vts(path):
    """Return (points (n, 3), pointdata {name: (n,)}, dims (nx, ny, nz))."""
    r = vtkXMLStructuredGridReader()
    r.SetFileName(path); r.Update()
    g = r.GetOutput()
    dims = [0, 0, 0]
    g.GetDimensions(dims)
    nx, ny, nz = dims
    pts = vtk_to_numpy(g.GetPoints().GetData()).copy()
    pd = g.GetPointData()
    fields = {}
    for k in range(pd.GetNumberOfArrays()):
        a = pd.GetArray(k)
        fields[a.GetName()] = vtk_to_numpy(a).copy()
    return pts, fields, (nx, ny, nz)


def write_vts(path, pts_3d, fields, dims):
    """Write a structured grid .vts. pts_3d: (n, 3); fields: {name: (n,)};
    dims: (nx, ny, nz). VTK uses Fortran ordering (i fastest, then j, then k)."""
    g = vtkStructuredGrid()
    g.SetDimensions(dims)
    pv = vtkPoints()
    pv.SetData(numpy_to_vtk(pts_3d.astype(np.float32), deep=1))
    g.SetPoints(pv)
    pd = g.GetPointData()
    for name, vals in fields.items():
        a = numpy_to_vtk(vals.astype(np.float32), deep=1)
        a.SetName(name)
        pd.AddArray(a)
    w = vtkXMLStructuredGridWriter()
    w.SetFileName(path); w.SetInputData(g)
    w.SetDataModeToAscii()
    w.Write()


def transform_section(pts_2d, dmp, delta_theta, x_spl, r_spl):
    """pts_2d: (n, 3) with cols (m', theta, 0). Returns (n, 3) cartesian."""
    m_blade = pts_2d[:, 0]; theta = pts_2d[:, 1]
    m_3D = m_blade + dmp                                          # delmp folded into dmp by caller
    x = np.asarray(x_spl(m_3D), dtype=np.float64)
    r = np.asarray(r_spl(m_3D), dtype=np.float64)
    th = theta + delta_theta
    out = np.column_stack([x, r * np.sin(th), r * np.cos(th)])
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument("sections_root", help="Dir containing section_1..N subdirs")
    p.add_argument("input_path", help="Path to 3dbgbinput.X.dat")
    p.add_argument("output_dir", help="Output dir for stacked .vts/.vtm")
    p.add_argument("--vts-name", default="final_grid.vts")
    p.add_argument("--le-idx", type=int, default=120,
                   help="Index of LE in the blade.X.X.SECTION file (i_LE_topo "
                        "from build_9block_topology). Default 120 for 240-pt blades.")
    p.add_argument("--blade-glob", default="blade.{ia}.*",
                   help="Glob inside section_{ia} for the blade file. Use {ia}.")
    args = p.parse_args()

    cfg = parse_3dbgbinput(args.input_path)
    nsl = cfg["nsl"]
    print(f"Parsed 3dbgbinput: nsl={nsl}, sweep_cps={cfg['sweep_cps'].shape}, "
          f"lean_cps={cfg['lean_cps'].shape}, "
          f"LE line {cfg['le_hub']}..{cfg['le_cas']}")

    # Per-section transformation parameters
    sec_data = []
    for ia in range(1, nsl + 1):
        sl = cfg["streams"][ia - 1]
        mp_s, x_spl, r_spl = streamline_splines(sl)
        msle = line_intersection_msle(mp_s, x_spl, r_spl,
                                       cfg["le_hub"], cfg["le_cas"])
        # span maps section index to [0, 1]
        span = (ia - 1) / max(1, nsl - 1)
        delmp = interp_cps(span, cfg["sweep_cps"])
        delta_theta = interp_cps(span, cfg["lean_cps"])

        # dmp = msle - m'_blade(LE_idx, ia). Read blade file to get LE m'.
        import glob
        sec_dir = os.path.join(args.sections_root, f"section_{ia}")
        blade_glob = args.blade_glob.format(ia=ia)
        blade_files = glob.glob(os.path.join(sec_dir, blade_glob))
        if not blade_files:
            raise FileNotFoundError(
                f"No blade file matching {blade_glob} in {sec_dir}")
        blade_raw = np.loadtxt(blade_files[0], skiprows=2)
        m_blade_LE = blade_raw[args.le_idx, 0]
        dmp = msle - m_blade_LE + delmp

        sec_data.append({
            "ia": ia, "x_spl": x_spl, "r_spl": r_spl, "msle": msle,
            "dmp": dmp, "delta_theta": delta_theta, "span": span,
        })
        if ia == 1 or ia == nsl:
            print(f"  section {ia}: msle={msle:.4f}, m'_blade_LE={m_blade_LE:.4f}, "
                  f"dmp={dmp:.4f}, delta_theta={delta_theta:.4e}")

    # For each block, stack 21 sections into a 3D grid
    os.makedirs(args.output_dir, exist_ok=True)
    block_files = []
    for label, subdir in BLOCKS:
        # Read section_1's block to determine (nx, ny)
        first = os.path.join(args.sections_root, "section_1", subdir, args.vts_name)
        pts0, fields0, dims0 = read_vts(first)
        nx, ny, _ = dims0
        n_xy = nx * ny

        all_pts = np.zeros((n_xy * nsl, 3), dtype=np.float64)
        all_fields = {name: np.zeros(n_xy * nsl, dtype=np.float64)
                       for name in fields0}

        for ia in range(1, nsl + 1):
            path = os.path.join(args.sections_root, f"section_{ia}", subdir,
                                 args.vts_name)
            pts_2d, fields, dims = read_vts(path)
            if dims != dims0:
                raise RuntimeError(
                    f"{label}: section_{ia} has dims {dims}, "
                    f"expected {dims0} (use unified mode for stackable grids).")
            sd = sec_data[ia - 1]
            cart = transform_section(pts_2d, sd["dmp"], sd["delta_theta"],
                                      sd["x_spl"], sd["r_spl"])
            offset = (ia - 1) * n_xy
            all_pts[offset:offset + n_xy] = cart
            for name, vals in fields.items():
                all_fields[name][offset:offset + n_xy] = vals

        out_name = f"{label.replace('.', '_')}.vts"
        out_path = os.path.join(args.output_dir, out_name)
        write_vts(out_path, all_pts, all_fields, (nx, ny, nsl))
        block_files.append((label, out_name))
        print(f"  wrote {out_path}: {nx}x{ny}x{nsl}")

    # VTM index
    blocks_xml = []
    for k, (label, fn) in enumerate(block_files):
        blocks_xml.append(
            f'    <Block index="{k}" name="{label}">\n'
            f'      <DataSet index="0" file="{fn}"/>\n'
            f'    </Block>')
    vtm = (f'<?xml version="1.0"?>\n'
           f'<VTKFile type="vtkMultiBlockDataSet" version="1.0" '
           f'byte_order="LittleEndian">\n'
           f'  <vtkMultiBlockDataSet>\n'
           + "\n".join(blocks_xml) + "\n"
           f'  </vtkMultiBlockDataSet>\n'
           f'</VTKFile>\n')
    vtm_path = os.path.join(args.output_dir, "stacked_3d.vtm")
    with open(vtm_path, "w") as f:
        f.write(vtm)
    print(f"\nwrote {vtm_path}")


if __name__ == "__main__":
    main()
