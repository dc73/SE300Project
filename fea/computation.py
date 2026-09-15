"""Read the project input format and produce demonstration result fields."""

from pathlib import Path


RESULT_COLUMNS = {
    "Stress X": (1, 0.2), "Stress Y": (2, 0.2),
    "Strain X": (1, 0.3), "Strain Y": (2, 0.3),
    "Von Mises Stress": (1, 0.4), "Integration": (1, 0.5),
}


def _data_lines(path):
    with Path(path).open(encoding="utf-8") as input_file:
        return [line.split("#", 1)[0].strip() for line in input_file
                if line.split("#", 1)[0].strip()]


def run_analysis(input_file="INPUT_FEA_PROTUS_3.txt"):
    """Read an FEA input file and return coordinates, elements, and demo results."""
    print("Processing FEA Analysis...")
    try:
        lines = _data_lines(input_file)
        if len(lines) < 11:
            raise ValueError("header is incomplete")
        num_nodes, num_elements, num_bcs, num_loads = map(int, lines[7:11])
        expected = 11 + num_nodes + num_elements + num_bcs + num_loads
        if len(lines) < expected:
            raise ValueError(f"expected {expected} data lines, found {len(lines)}")
        element_start = 11 + num_nodes
        nodes = [[float(value) for value in row.split(",")]
                 for row in lines[11:element_start]]
        elements = [[int(value) for value in row.split(",")[:5]]
                    for row in lines[element_start:element_start + num_elements]]
        if any(len(node) != 3 for node in nodes):
            raise ValueError("each node must contain id,x,y")
        if any(len(element) != 5 for element in elements):
            raise ValueError("each element must contain id and four node ids")
    except (OSError, ValueError) as error:
        print(f"Could not read '{input_file}': {error}")
        return None, None, None

    node_coords = [(node[1], node[2]) for node in nodes]
    # Placeholder fields for visualization; this is not yet a physical FEA solver.
    results = {name: [node[index] * scale for node in nodes]
               for name, (index, scale) in RESULT_COLUMNS.items()}
    print(f"Loaded {num_nodes} nodes and {num_elements} elements.")
    return node_coords, elements, results
