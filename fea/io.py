"""Persistence helpers for FEA conversion and result files."""

import json
from pathlib import Path


def write_conversion_file(node_list, element_list, bc_list, load_list,
                          material, analysis_type,
                          filename="INPUT_FEA_PROTUS_3.txt"):
    sections = [
        ("MATERIAL DATA", [analysis_type, material["E"], material["Poisson"],
                           material["Density"], 1.0, 0.0, -9.81, len(node_list),
                           len(element_list), len(bc_list), len(load_list)]),
        ("NODE DATA", node_list), ("ELEMENT DATA", element_list),
        ("BC, 1-Fixed or 0-Free", bc_list), ("LOAD", load_list),
    ]
    try:
        with Path(filename).open("w", encoding="utf-8") as output:
            for title, rows in sections:
                output.write(f"#---------------------{title}\n")
                for row in rows:
                    output.write(",".join(map(str, row)) + "\n"
                                 if isinstance(row, (list, tuple)) else f"{row}\n")
        return True
    except OSError as error:
        print(f"Error writing conversion file: {error}")
        return False


def export_file(data, save_path):
    try:
        with Path(save_path).open("w", encoding="utf-8") as output:
            json.dump(data, output, indent=2)
        return True
    except (OSError, TypeError) as error:
        print(f"Error exporting file: {error}")
        return False


def load_file(load_path):
    try:
        with Path(load_path).open(encoding="utf-8") as input_file:
            return json.load(input_file)
    except (OSError, json.JSONDecodeError) as error:
        print(f"Error loading file: {error}")
        return None
