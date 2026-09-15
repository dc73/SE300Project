"""Convert planar mesh geometry to the project's node and element records."""


def build_mesh_data(vertices, faces):
    """Return nodes, supported faces, and nodes fixed along the left edge."""
    nodes = [(index, float(point[0]), float(point[1]))
             for index, point in enumerate(vertices, start=1)]
    if not nodes:
        raise ValueError("The selected file does not contain any vertices.")

    elements = []
    for face in faces:
        node_ids = [int(index) + 1 for index in face]
        if len(node_ids) == 3:
            node_ids.append(node_ids[-1])
        if len(node_ids) == 4:
            elements.append((len(elements) + 1, *node_ids))
    if not elements:
        raise ValueError("The mesh contains no triangular or quadrilateral faces.")

    min_x = min(node[1] for node in nodes)
    span = max(node[1] for node in nodes) - min_x
    tolerance = max(span * 1e-9, 1e-12)
    fixed_nodes = [node[0] for node in nodes if abs(node[1] - min_x) <= tolerance]
    return nodes, elements, [(node_id, 1, 1) for node_id in fixed_nodes]
