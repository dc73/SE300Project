import re
import sys

class FEAComputation:
    """Handles FEA computations and matrix operations."""
    
    @staticmethod
    def matrix_multiply(A, B):
        RowsA = len(A)
        ColsA = len(A[0])
        ColsB = len(B[0])
        Result = [[0 for _ in range(ColsB)] for _ in range(RowsA)]
        for i in range(RowsA):
            for j in range(ColsB):
                for k in range(ColsA):
                    Result[i][j] += A[i][k] * B[k][j]
        return Result
    
    @staticmethod
    def matrix_add(A, B):
        Rows = len(A)
        Cols = len(A[0])
        Result = [[0 for _ in range(Cols)] for _ in range(Rows)]
        for i in range(Rows):
            for j in range(Cols):
                Result[i][j] = A[i][j] + B[i][j]
        return Result

    @staticmethod
    def matrix_scalar_multiply(A, scalar):
        Rows = len(A)
        Cols = len(A[0])
        return [[A[i][j] * scalar for j in range(Cols)] for i in range(Rows)]
    
    @staticmethod
    def matrix_transpose(A):
        Rows = len(A)
        Cols = len(A[0])
        return [[A[i][j] for i in range(Rows)] for j in range(Cols)]
    
    @staticmethod
    def run_analysis(input_file="INPUT_FEA_PROTUS_3.txt"):
        """
        Reads the FEA input file, performs dummy analysis,
        and returns node coordinates, element data, and a dictionary of results.
        """
        print("*******************************************************************************")
        print("Processing FEA Analysis...")
        try:
            with open(input_file, "r") as f:
                InputLines = f.readlines()
        except FileNotFoundError:
            print(f"Error: '{input_file}' not found. Please convert a mesh first.")
            return None, None, None
        
        # Read header parameters
        AnalysisType = int(re.split(r'\s+', InputLines[1])[0])
        ElasticModulus = float(re.split(r'\s+', InputLines[2])[0])
        PoissonsRatio = float(re.split(r'\s+', InputLines[3])[0])
        MaterialDensity = float(re.split(r'\s+', InputLines[4])[0])
        Thickness = float(re.split(r'\s+', InputLines[5])[0])
        AccelX = float(re.split(r'\s+', InputLines[6])[0])
        AccelY = float(re.split(r'\s+', InputLines[7])[0])
        NumNodes = int(re.split(r'\s+', InputLines[8])[0])
        NumElements = int(re.split(r'\s+', InputLines[9])[0])
        NumBC = int(re.split(r'\s+', InputLines[10])[0])
        NumLoads = int(re.split(r'\s+', InputLines[11])[0])
        
        DOFTotal = 2 * NumNodes
        DOFPerNode = 2
        
        print("END OF USER INPUT PROCESSING")
        print("================= SOLVER ======================")
        print("NUMBER OF NODES IS -", NumNodes)
        print("NUMBER OF ELEMENTS IS -", NumElements)
        print("TOTAL NUMBER OF VARIABLES IN THE MODEL -", DOFTotal)
        
        # Read Node Data
        NodeData = [[0, 0, 0] for _ in range(NumNodes)]
        Index = 13
        while Index < (NumNodes + 13):
            NodeValues = [float(val.strip()) for val in re.split(r',', InputLines[Index])]
            NodeData[Index - 13][0] = NodeValues[0]  # node id
            NodeData[Index - 13][1] = NodeValues[1]  # x coordinate
            NodeData[Index - 13][2] = NodeValues[2]  # y coordinate
            Index += 1
        
        # Read Element Data
        ElementData = []
        Index = 14 + NumNodes
        while Index < (14 + NumNodes + NumElements):
            ElementValues = [int(val.strip()) for val in re.split(r',', InputLines[Index])]
            ElementData.append(ElementValues[:5])
            Index += 1
        
        # Dummy Post-Processing: Generate result values for demonstration
        NumNodes = len(NodeData)
        ResultsMatrix = [[None] * 72 for _ in range(NumNodes)]
        for i in range(NumNodes):
            ResultsMatrix[i][1] = NodeData[i][1] * 0.2  # Stress X
            ResultsMatrix[i][6] = NodeData[i][2] * 0.2  # Stress Y
            ResultsMatrix[i][66] = NodeData[i][1] * 0.3  # Strain X
            ResultsMatrix[i][70] = NodeData[i][2] * 0.3  # Strain Y
            ResultsMatrix[i][30] = NodeData[i][1] * 0.4  # Von Mises Stress
            ResultsMatrix[i][15] = NodeData[i][1] * 0.5  # Integration dummy
        
        results = {
            "Stress X": [ResultsMatrix[i][1] for i in range(NumNodes)],
            "Stress Y": [ResultsMatrix[i][6] for i in range(NumNodes)],
            "Strain X": [ResultsMatrix[i][66] for i in range(NumNodes)],
            "Strain Y": [ResultsMatrix[i][70] for i in range(NumNodes)],
            "Von Mises Stress": [ResultsMatrix[i][30] for i in range(NumNodes)],
            "Integration": [ResultsMatrix[i][15] for i in range(NumNodes)]
        }
        
        print("FINISHED")
        node_coords = [(NodeData[i][1], NodeData[i][2]) for i in range(NumNodes)]
        return node_coords, ElementData, results
