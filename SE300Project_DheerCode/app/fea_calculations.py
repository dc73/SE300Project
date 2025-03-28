import tkinter as tk
from tkinter import messagebox
import numpy as np
import re
import math
import io
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.tri as mtri

# ------------------------- MATRIX OPERATIONS -------------------------
def MatrixMultiply(MatrixA, MatrixB):
    RowsA = len(MatrixA)
    ColsA = len(MatrixA[0])
    ColsB = len(MatrixB[0])
    Result = [[0 for _ in range(ColsB)] for _ in range(RowsA)]
    for i in range(RowsA):
        for j in range(ColsB):
            for k in range(ColsA):
                Result[i][j] += MatrixA[i][k] * MatrixB[k][j]
    return Result

def MatrixAdd(MatrixA, MatrixB):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    Result = [[0 for _ in range(Cols)] for _ in range(Rows)]
    for i in range(Rows):
        for j in range(Cols):
            Result[i][j] = MatrixA[i][j] + MatrixB[i][j]
    return Result

def MatrixScalarMultiply(MatrixA, Scalar):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    return [[MatrixA[i][j] * Scalar for j in range(Cols)] for i in range(Rows)]

def MatrixTranspose(MatrixA):
    Rows = len(MatrixA)
    Cols = len(MatrixA[0])
    return [[MatrixA[i][j] for i in range(Rows)] for j in range(Cols)]

# ----- Local helper functions (for element-level computations) -----
def JacobianDeterminant_local(XCoords, YCoords, s, t):
    J00 = (t-1)*XCoords[0][0] + (1-t)*XCoords[0][1] + (1+t)*XCoords[0][2] - (1+t)*XCoords[0][3]
    J01 = (t-1)*YCoords[0][0] + (1-t)*YCoords[0][1] + (1+t)*YCoords[0][2] - (1+t)*YCoords[0][3]
    J10 = (s-1)*XCoords[0][0] - (1+s)*XCoords[0][1] + (1+s)*XCoords[0][2] + (1-s)*XCoords[0][3]
    J11 = (s-1)*YCoords[0][0] - (1+s)*YCoords[0][1] + (1+s)*YCoords[0][2] + (1-s)*YCoords[0][3]
    return (J00 * J11 - J01 * J10) / 16.0

def ComputeBMatrix_local(XCoords, YCoords, s, t):
    B = [[0]*8 for _ in range(3)]
    a = 0.25 * (YCoords[0][0]*(s-1) + YCoords[0][1]*(-1-s) + YCoords[0][2]*(1+s) + YCoords[0][3]*(1-s))
    b = 0.25 * (YCoords[0][0]*(t-1) + YCoords[0][1]*(1-t) + YCoords[0][2]*(1+t) + YCoords[0][3]*(-1-t))
    c = 0.25 * (XCoords[0][0]*(t-1) + XCoords[0][1]*(1-t) + XCoords[0][2]*(1+t) + XCoords[0][3]*(-1-t))
    d = 0.25 * (XCoords[0][0]*(s-1) + XCoords[0][1]*(-1-s) + XCoords[0][2]*(1+s) + XCoords[0][3]*(1-s))
    Jdet = JacobianDeterminant_local(XCoords, YCoords, s, t)
    B[0][0] = (a*0.25*(t-1) - b*0.25*(s-1)) / Jdet
    B[2][0] = (c*0.25*(s-1) - d*0.25*(t-1)) / Jdet
    B[1][1] = (c*0.25*(s-1) - d*0.25*(t-1)) / Jdet
    B[2][1] = (a*0.25*(t-1) - b*0.25*(s-1)) / Jdet
    B[0][2] = (a*0.25*(1-t) - b*0.25*(-1-s)) / Jdet
    B[2][2] = (c*0.25*(-1-s) - d*0.25*(1-t)) / Jdet
    B[1][3] = (c*0.25*(-1-s) - d*0.25*(1-t)) / Jdet
    B[2][3] = (a*0.25*(1-t) - b*0.25*(-1-s)) / Jdet
    B[0][4] = (a*0.25*(1+t) - b*0.25*(1+s)) / Jdet
    B[2][4] = (c*0.25*(1+s) - d*0.25*(1+t)) / Jdet
    B[1][5] = (c*0.25*(1+s) - d*0.25*(1+t)) / Jdet
    B[2][5] = (a*0.25*(1+t) - b*0.25*(1+s)) / Jdet
    B[0][6] = (a*0.25*(-1-t) - b*0.25*(1-s)) / Jdet
    B[2][6] = (c*0.25*(1-s) - d*0.25*(-1-t)) / Jdet
    B[1][7] = (c*0.25*(1-s) - d*0.25*(-1-t)) / Jdet
    B[2][7] = (a*0.25*(-1-t) - b*0.25*(1-s)) / Jdet
    return B

# ------------------------- FEA ANALYSIS FUNCTION -------------------------
def RunAnalysis():
    """
    Reads the FEA input file "INPUT_FEA.txt", performs FEA analysis,
    assembles global stiffness and mass matrices, applies boundary conditions,
    solves for nodal displacements, performs dummy post-processing, and
    writes a VTK file "SE300Paraview.vtk" for Paraview.
    
    Returns:
      - node_coords: list of (x,y) nodal coordinates,
      - ElementData: connectivity data,
      - results: a dictionary with keys:
            "Nodal Displacements", "Stress X", "Stress Y", "Strain X", "Strain Y", "Von Mises Stress"
    """
    print("*******************************************************************************")
    print("Processing FEA Analysis...")
    
    try:
        InputFile = open('INPUT_FEA.txt', 'r')
    except FileNotFoundError:
        print("Error: 'INPUT_FEA.txt' not found. Please convert a mesh first.")
        return None, None, None
    InputLines = InputFile.readlines()
    InputFile.close()
    
    AnalysisType    = int(re.split('\s+', InputLines[1])[0])
    ElasticModulus  = float(re.split('\s+', InputLines[2])[0])
    PoissonsRatio   = float(re.split('\s+', InputLines[3])[0])
    MaterialDensity = float(re.split('\s+', InputLines[4])[0])
    Thickness       = float(re.split('\s+', InputLines[5])[0])
    AccelX          = float(re.split('\s+', InputLines[6])[0])
    AccelY          = float(re.split('\s+', InputLines[7])[0])
    NumNodes        = int(re.split('\s+', InputLines[8])[0])
    NumElements     = int(re.split('\s+', InputLines[9])[0])
    NumBC           = int(re.split('\s+', InputLines[10])[0])
    NumLoads        = int(re.split('\s+', InputLines[11])[0])
    
    DOFTotal        = 2 * NumNodes
    DOFPerNode      = 2
    Pi              = np.arccos(-1)
    
    if AnalysisType == 1:
        print("Type of Analysis: Plane Stress")
    else:
        print("Type of Analysis: Plane Strain")
    
    NodeData = [[0 for _ in range(3)] for _ in range(NumNodes)]
    Index = 13
    while Index < (NumNodes + 13):
        NodeValues = [float(val.strip()) for val in re.split(',', InputLines[Index])]
        NodeData[Index - 13] = NodeValues[:3]
        Index += 1
    
    ElementData = [[0 for _ in range(5)] for _ in range(NumElements)]
    Index = 14 + NumNodes
    while Index < (14 + NumNodes + NumElements):
        ElementValues = [int(val.strip()) for val in re.split(',', InputLines[Index])]
        ElementData[Index - (14 + NumNodes)] = ElementValues[:5]
        Index += 1
    
    BCVector = [[0] for _ in range(DOFTotal)]
    Index = 15 + NumNodes + NumElements
    while Index < (15 + NumNodes + NumElements + NumBC):
        BCValues = [int(val.strip()) for val in re.split(',', InputLines[Index])]
        BCVector[2 * BCValues[0] - 2][0] = BCValues[1]
        BCVector[2 * BCValues[0] - 1][0] = BCValues[2]
        Index += 1
    
    LoadsVector = [[0] for _ in range(DOFTotal)]
    Index = 16 + NumNodes + NumElements + NumBC
    while Index < (16 + NumNodes + NumElements + NumBC + NumLoads):
        LoadValues = [val.strip() for val in re.split(',', InputLines[Index])]
        LoadsVector[2 * int(LoadValues[0]) - 2][0] = float(LoadValues[1])
        LoadsVector[2 * int(LoadValues[0]) - 1][0] = float(LoadValues[2])
        Index += 1
    
    AccelerationVector = [[0] for _ in range(DOFTotal)]
    for r in range(0, DOFTotal, 2):
        AccelerationVector[r][0] = AccelX
    for r in range(1, DOFTotal, 2):
        AccelerationVector[r][0] = AccelY
    
    print("END OF USER INPUT PROCESSING")
    print("================= SOLVER ======================")
    print("NUMBER OF NODES IS -", NumNodes)
    print("NUMBER OF ELEMENTS IS -", NumElements)
    print("TOTAL NUMBER OF VARIABLES IN THE MODEL -", DOFTotal)
    
    s1, t1 = -0.577350269, -0.577350269
    s2, t2 =  0.577350269, -0.577350269
    s3, t3 =  0.577350269,  0.577350269
    s4, t4 = -0.577350269,  0.577350269
    
    GlobalStiffness = [[0 for _ in range(DOFTotal)] for _ in range(DOFTotal)]
    GlobalMass      = [[0 for _ in range(DOFTotal)] for _ in range(DOFTotal)]
    
    
    ParaviewResults = [[None for _ in range(72)] for _ in range(NumNodes)]
    for i in range(NumNodes):
        ParaviewResults[i][0] = NodeData[i][0]
    
    if AnalysisType == 1:
        DMatrix = [[ElasticModulus/(1-PoissonsRatio**2), (PoissonsRatio*ElasticModulus)/(1-PoissonsRatio**2), 0],
                   [(PoissonsRatio*ElasticModulus)/(1-PoissonsRatio**2), ElasticModulus/(1-PoissonsRatio**2), 0],
                   [0, 0, (ElasticModulus*(1-PoissonsRatio))/(2*(1-PoissonsRatio**2))]]
    else:
        DMatrix = [[ElasticModulus/((1+PoissonsRatio)*(1-2*PoissonsRatio))*(1-PoissonsRatio),
                    PoissonsRatio*ElasticModulus/((1+PoissonsRatio)*(1-2*PoissonsRatio)), 0],
                   [PoissonsRatio*ElasticModulus/((1+PoissonsRatio)*(1-2*PoissonsRatio)),
                    ElasticModulus/((1+PoissonsRatio)*(1-2*PoissonsRatio))*(1-PoissonsRatio), 0],
                   [0, 0, ElasticModulus/((1+PoissonsRatio)*(1-2*PoissonsRatio))*((1-2*PoissonsRatio)/2)]]
    
    for elem in range(1, NumElements+1):
        XCoords = [[ NodeData[ ElementData[elem-1][1]-1 ][1],
                     NodeData[ ElementData[elem-1][2]-1 ][1],
                     NodeData[ ElementData[elem-1][3]-1 ][1],
                     NodeData[ ElementData[elem-1][4]-1 ][1] ]]
        YCoords = [[ NodeData[ ElementData[elem-1][1]-1 ][2],
                     NodeData[ ElementData[elem-1][2]-1 ][2],
                     NodeData[ ElementData[elem-1][3]-1 ][2],
                     NodeData[ ElementData[elem-1][4]-1 ][2] ]]
        
        Ke1 = MatrixScalarMultiply(
                MatrixScalarMultiply(
                  MatrixMultiply(
                    MatrixTranspose(ComputeBMatrix_local(XCoords, YCoords, s1, t1)),
                    MatrixMultiply(DMatrix, ComputeBMatrix_local(XCoords, YCoords, s1, t1))
                  ),
                  JacobianDeterminant_local(XCoords, YCoords, s1, t1)
                ),
                Thickness)
        Ke2 = MatrixScalarMultiply(
                MatrixScalarMultiply(
                  MatrixMultiply(
                    MatrixTranspose(ComputeBMatrix_local(XCoords, YCoords, s2, t2)),
                    MatrixMultiply(DMatrix, ComputeBMatrix_local(XCoords, YCoords, s2, t2))
                  ),
                  JacobianDeterminant_local(XCoords, YCoords, s2, t2)
                ),
                Thickness)
        Ke3 = MatrixScalarMultiply(
                MatrixScalarMultiply(
                  MatrixMultiply(
                    MatrixTranspose(ComputeBMatrix_local(XCoords, YCoords, s3, t3)),
                    MatrixMultiply(DMatrix, ComputeBMatrix_local(XCoords, YCoords, s3, t3))
                  ),
                  JacobianDeterminant_local(XCoords, YCoords, s3, t3)
                ),
                Thickness)
        Ke4 = MatrixScalarMultiply(
                MatrixScalarMultiply(
                  MatrixMultiply(
                    MatrixTranspose(ComputeBMatrix_local(XCoords, YCoords, s4, t4)),
                    MatrixMultiply(DMatrix, ComputeBMatrix_local(XCoords, YCoords, s4, t4))
                  ),
                  JacobianDeterminant_local(XCoords, YCoords, s4, t4)
                ),
                Thickness)
        Ke = MatrixAdd(MatrixAdd(MatrixAdd(Ke1, Ke2), Ke3), Ke4)
        
        def ComputeNMatrix(s, t):
            NMatrix = [[0]*8, [0]*8]
            NMatrix[0][0] = ((1-s)*(1-t)/4)
            NMatrix[1][1] = ((1-s)*(1-t)/4)
            NMatrix[0][2] = ((1+s)*(1-t)/4)
            NMatrix[1][3] = ((1+s)*(1-t)/4)
            NMatrix[0][4] = ((1+s)*(1+t)/4)
            NMatrix[1][5] = ((1+s)*(1+t)/4)
            NMatrix[0][6] = ((1-s)*(1+t)/4)
            NMatrix[1][7] = ((1-s)*(1+t)/4)
            return NMatrix
        
        MassMatrix1 = MatrixScalarMultiply(
                        MatrixScalarMultiply(
                          MatrixScalarMultiply(
                            MatrixMultiply(MatrixTranspose(ComputeNMatrix(s1, t1)), ComputeNMatrix(s1, t1)),
                            JacobianDeterminant_local(XCoords, YCoords, s1, t1)
                          ),
                          Thickness
                        ),
                        MaterialDensity)
        MassMatrix2 = MatrixScalarMultiply(
                        MatrixScalarMultiply(
                          MatrixScalarMultiply(
                            MatrixMultiply(MatrixTranspose(ComputeNMatrix(s2, t2)), ComputeNMatrix(s2, t2)),
                            JacobianDeterminant_local(XCoords, YCoords, s2, t2)
                          ),
                          Thickness
                        ),
                        MaterialDensity)
        MassMatrix3 = MatrixScalarMultiply(
                        MatrixScalarMultiply(
                          MatrixScalarMultiply(
                            MatrixMultiply(MatrixTranspose(ComputeNMatrix(s3, t3)), ComputeNMatrix(s3, t3)),
                            JacobianDeterminant_local(XCoords, YCoords, s3, t3)
                          ),
                          Thickness
                        ),
                        MaterialDensity)
        MassMatrix4 = MatrixScalarMultiply(
                        MatrixScalarMultiply(
                          MatrixScalarMultiply(
                            MatrixMultiply(MatrixTranspose(ComputeNMatrix(s4, t4)), ComputeNMatrix(s4, t4)),
                            JacobianDeterminant_local(XCoords, YCoords, s4, t4)
                          ),
                          Thickness
                        ),
                        MaterialDensity)
        MassMatrixElement = MatrixAdd(MatrixAdd(MatrixAdd(MassMatrix1, MassMatrix2), MassMatrix3), MassMatrix4)
        
        for i in range(0, 4):
            for j in range(0, 2):
                for k in range(0, 4):
                    for m in range(0, 2):
                        RowIndex = DOFPerNode * (ElementData[elem-1][i+1] - 1) + j
                        ColIndex = DOFPerNode * (ElementData[elem-1][k+1] - 1) + m
                        GlobalStiffness[RowIndex][ColIndex] += Ke[DOFPerNode*i+j][DOFPerNode*k+m]
                        GlobalMass[RowIndex][ColIndex] += MassMatrixElement[DOFPerNode*i+j][DOFPerNode*k+m]
    
    print("GLOBAL STIFFNESS MATRIX READY")
    
    CombinedLoads = np.add(MatrixMultiply(MatrixTranspose(GlobalMass), AccelerationVector), LoadsVector)
    for i in range(DOFTotal):
        if BCVector[i][0] == 1:
            CombinedLoads[i][0] = 0
    
    ModifiedStiffness = [[0 for _ in range(DOFTotal)] for _ in range(DOFTotal)]
    for i in range(DOFTotal):
        if BCVector[i][0] == 1:
            for j in range(DOFTotal):
                ModifiedStiffness[i][j] = 0
                ModifiedStiffness[j][i] = 0
            ModifiedStiffness[i][i] = 1
        else:
            for j in range(i, DOFTotal):
                ModifiedStiffness[i][j] = GlobalStiffness[i][j]
            for j in range(i, DOFTotal):
                ModifiedStiffness[j][i] = GlobalStiffness[j][i]
    print("BOUNDARY CONDITIONS APPLIED")
    print("INVERTING GLOBAL STIFFNESS MATRIX")
    
    StiffnessInverse = np.linalg.pinv(ModifiedStiffness)
    print("GLOBAL STIFFNESS MATRIX INVERTED")
    
    NodalDisplacements = MatrixMultiply(StiffnessInverse, CombinedLoads)
    print("SYSTEM OF EQUATIONS SOLUTION IS DONE!")
    print("================ POST-PROCESSOR =================")
    print("STARTING POST-PROCESS DATA")
    
    # Use node x-coordinate to simulate variation in dummy results
    nodal_result = [NodeData[i][1] for i in range(NumNodes)]
    
    results = {
        "Nodal Displacements": nodal_result,
        "Stress X": [val * 0.1 for val in nodal_result],
        "Stress Y": [val * 0.2 for val in nodal_result],
        "Strain X": [val * 0.3 for val in nodal_result],
        "Strain Y": [val * 0.4 for val in nodal_result],
        "Von Mises Stress": [val * 0.5 for val in nodal_result]
    }
    
    node_coords = [(NodeData[i][1], NodeData[i][2]) for i in range(NumNodes)]
    
    # Write VTK file for Paraview
    vtk_file = open('SE300Paraview.vtk', 'w')
    vtk_file.write('# vtk DataFile Version 2.0\n')
    vtk_file.write('FEA Results\n')
    vtk_file.write('ASCII\n')
    vtk_file.write('DATASET UNSTRUCTURED_GRID\n')
    vtk_file.write('POINTS ' + str(NumNodes) + ' float\n')
    for i in range(NumNodes):
        vtk_file.write(str(NodeData[i][1]) + ' ' + str(NodeData[i][2]) + ' 0\n')
    vtk_file.write('CELLS ' + str(NumElements) + ' ' + str(5*NumElements) + '\n')
    for i in range(NumElements):
        vtk_file.write('4 ' + ' '.join([str(ElementData[i][j]-1) for j in range(1,5)]) + '\n')
    vtk_file.write('CELL_TYPES ' + str(NumElements) + '\n')
    for i in range(NumElements):
        vtk_file.write('9\n')
    vtk_file.write('POINT_DATA ' + str(NumNodes) + '\n')
    vtk_file.write('SCALARS Nodal_Displacement float\n')
    vtk_file.write('LOOKUP_TABLE default\n')
    for i in range(NumNodes):
        vtk_file.write(str(nodal_result[i]) + '\n')
    vtk_file.close()
    print("VTK FILE 'SE300Paraview.vtk' CREATED FOR PARAVIEW")
    
    return node_coords, ElementData, results

# ------------------------- PLOTTING FUNCTIONS -------------------------
def plot_2d_result(parent, node_coords, element_data, result_values, result_label):
    # If result_values length equals 2*number of nodes, average them (not needed for our dummy data)
    if len(result_values) == 2 * len(node_coords):
        averaged = []
        for i in range(len(node_coords)):
            avg_val = (result_values[2*i] + result_values[2*i+1]) / 2
            averaged.append(avg_val)
        result_values = averaged

    fig, ax = plt.subplots(figsize=(6, 6))
    fig.patch.set_facecolor('#4A4A4A')
    ax.set_facecolor('#4A4A4A')
    x_vals = [pt[0] for pt in node_coords]
    y_vals = [pt[1] for pt in node_coords]
    triangles = []
    for elem in element_data:
        indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
        if indices[2] == indices[3]:
            triangles.append([indices[0], indices[1], indices[2]])
        else:
            triangles.append([indices[0], indices[1], indices[2]])
            triangles.append([indices[0], indices[2], indices[3]])
    triangulation = mtri.Triangulation(x_vals, y_vals, triangles)
    tpc = ax.tripcolor(triangulation, result_values, shading='gouraud', cmap='viridis')
    for elem in element_data:
        indices = [elem[1] - 1, elem[2] - 1, elem[3] - 1, elem[4] - 1]
        poly_x = [node_coords[i][0] for i in indices]
        poly_y = [node_coords[i][1] for i in indices]
        poly_x.append(poly_x[0])
        poly_y.append(poly_y[0])
        ax.plot(poly_x, poly_y, color='black', linewidth=1)
    ax.set_title(f"{result_label} (2D)", color='white')
    ax.set_xlabel("X Axis", color='white')
    ax.set_ylabel("Y Axis", color='white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.set_aspect('equal', 'box')
    cbar = fig.colorbar(tpc, ax=ax)
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.outline.set_edgecolor('white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.draw()
    widget = canvas.get_tk_widget()
    widget.pack(fill='both', expand=True)
    return canvas
