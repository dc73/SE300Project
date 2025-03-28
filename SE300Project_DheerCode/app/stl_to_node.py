import trimesh

def convert_stl_to_node(mesh_path, mesh_density, material_props, analysis_type, force_node, force_value):
    """
    Converts a mesh file into an FEA input file.
    
    :param mesh_path: Path to the input mesh (STL/OBJ) file.
    :param mesh_density: Integer (1 to 10) representing the mesh density (1=coarse, 10=fine).
    :param material_props: Dictionary with keys 'E', 'Poisson', and 'Density'.
    :param analysis_type: 1 for Plane Stress or 2 for Plane Strain.
    :param force_node: Node number where force is applied.
    :param force_value: Force value (N) at the specified node.
    :return: The name of the generated FEA input file.
    """
    mesh = trimesh.load(mesh_path)
    target_faces = max(100, int(len(mesh.faces) * (0.1 * mesh_density)))
    if hasattr(mesh, 'simplify_quadratic_decimation'):
        simplified = mesh.simplify_quadratic_decimation(target_faces)
    elif hasattr(mesh, 'simplify_vertex_clustering'):
        bbox = mesh.bounding_box.extents
        voxel_size = min(bbox) / (mesh_density * 5)
        simplified = mesh.simplify_vertex_clustering(voxel_size)
    else:
        simplified = mesh

    vertices = simplified.vertices
    nodes2d = vertices[:, :2]
    node_list = [(i + 1, coord[0], coord[1]) for i, coord in enumerate(nodes2d)]
    
    faces = simplified.faces
    element_list = []
    for i, face in enumerate(faces):
        face = list(face)
        if len(face) == 3:
            element_list.append((i + 1, face[0] + 1, face[1] + 1, face[2] + 1, face[2] + 1))
        elif len(face) == 4:
            element_list.append((i + 1, face[0] + 1, face[1] + 1, face[2] + 1, face[3] + 1))
    
    max_node_id = max([node[0] for node in node_list])
    bc_list = [(nid, 1, 1) for nid in range(1, min(14, max_node_id + 1))]
    
    load_list = []
    try:
        force_node = int(force_node)
        force_value = float(force_value)
        if 1 <= force_node <= max_node_id:
            load_list.append((force_node, force_value, 1))
    except Exception:
        pass

    output_filename = "INPUT_FEA.txt"
    with open(output_filename, "w") as f:
        f.write("#---------------------MATERIAL DATA\n")
        f.write(f"{analysis_type}\n")
        f.write(f"{material_props['E']}\n")
        f.write(f"{material_props['Poisson']}\n")
        f.write(f"{material_props['Density']}\n")
        f.write("1.0\n")
        f.write("0.0\n")
        f.write("-9.81\n")
        f.write(f"{len(node_list)}\n")
        f.write(f"{len(element_list)}\n")
        f.write(f"{len(bc_list)}\n")
        f.write(f"{len(load_list)}\n")
        f.write("#---------------------NODE DATA\n")
        for node in node_list:
            f.write(f"{node[0]},{node[1]},{node[2]}\n")
        f.write("#---------------------ELEMENT DATA\n")
        for elem in element_list:
            f.write(f"{elem[0]},{elem[1]},{elem[2]},{elem[3]},{elem[4]}\n")
        f.write("#---------------------BC, 1-Fixed or 0-Free\n")
        for bc in bc_list:
            f.write(f"{bc[0]},{bc[1]},{bc[2]}\n")
        f.write("#---------------------LOAD AND BC\n")
        for load in load_list:
            f.write(f"{load[0]},{load[1]},{load[2]}\n")
    return output_filename
