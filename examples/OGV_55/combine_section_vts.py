"""Combine the 9 final_grid.vts files of one blade section into a single VTM
multi-block dataset. Preserves grid quality fields (orthogonality, aspect_ratio,
jacobian) since they live inside the .vts files unchanged.

Usage:
    python3 combine_section_vts.py <section_dir>

Writes <section_dir>/<section_name>.vtm referencing the 9 .vts files using
relative paths so the .vtm + .vts directory tree can be moved together.
"""
import os
import sys
import argparse


BLOCKS = [
    ("ogrid_3",                "ogrid_3"),
    ("top.left",               "top_passage_3/left"),
    ("top.center",             "top_passage_3/center"),
    ("top.right",              "top_passage_3/right"),
    ("bot.left",               "bot_passage_3/left"),
    ("bot.center",             "bot_passage_3/center"),
    ("bot.right",              "bot_passage_3/right"),
    ("middle.left",            "middle_passage_3/left"),
    ("middle.right",           "middle_passage_3/right"),
]


def write_vtm(section_dir, vts_name="final_grid.vts"):
    section_dir = os.path.abspath(section_dir)
    section_name = os.path.basename(section_dir.rstrip(os.sep))
    out_path = os.path.join(section_dir, f"{section_name}.vtm")

    missing = []
    blocks_xml = []
    for idx, (label, subdir) in enumerate(BLOCKS):
        rel = os.path.join(subdir, vts_name)
        full = os.path.join(section_dir, rel)
        if not os.path.isfile(full):
            missing.append(rel)
            continue
        blocks_xml.append(
            f'    <Block index="{idx}" name="{label}">\n'
            f'      <DataSet index="0" file="{rel}"/>\n'
            f'    </Block>')

    if missing:
        raise FileNotFoundError(
            f"section {section_name}: missing .vts files:\n  "
            + "\n  ".join(missing))

    body = "\n".join(blocks_xml)
    xml = (f'<?xml version="1.0"?>\n'
           f'<VTKFile type="vtkMultiBlockDataSet" version="1.0" '
           f'byte_order="LittleEndian">\n'
           f'  <vtkMultiBlockDataSet>\n'
           f'{body}\n'
           f'  </vtkMultiBlockDataSet>\n'
           f'</VTKFile>\n')

    with open(out_path, "w") as f:
        f.write(xml)
    return out_path


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument("section_dir", help="Directory containing the 9 block subdirs")
    p.add_argument("--vts-name", default="final_grid.vts",
                   help="VTS filename inside each block dir (default: final_grid.vts)")
    args = p.parse_args()

    out = write_vtm(args.section_dir, vts_name=args.vts_name)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
