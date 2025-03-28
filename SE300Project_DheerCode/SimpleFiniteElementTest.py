import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

nu = float(input("Input the possion's ratio: "))
E = float(input("Input modulus of elasticity [GPa]: "))

RectHeight = float(input("Input the height of the rectangle structure [m]: "))
RectLength = float(input("Input the length of the rectangle structure [m]: "))
Force = float(input("Input single force's magnitude in [N]: "))
MaterialThickness = float(input("Input the thickness of the material [m]: "))

ForceXComp = float(input("The X coordinate of the force on the beam [m]: "))
while (ForceXComp < 0 or ForceXComp > RectLength):
    print(f"Invalid value! Please enter a value between 0 and {RectLength:.3f}.")
    ForceXComp = float(input("The X coordinate of the force on the beam [m]: "))

print("The X component of the force is logged.")

# -------------------------------
# Calculation Section and Stiffness Matrix Calculation
# -------------------------------
Area1 = 0.5 * RectHeight * RectLength

# Element 1
LocalNodeX1E1 = 0.0
LocalNodeY1E1 = 0.0
LocalNodeX2E1 = RectLength
LocalNodeY2E1 = RectHeight
LocalNodeX3E1 = 0.0
LocalNodeY3E1 = RectHeight

b1E1 = LocalNodeY2E1 - LocalNodeY3E1
b2E1 = LocalNodeY3E1 - LocalNodeY1E1
b3E1 = LocalNodeY1E1 - LocalNodeY2E1

g1E1 = LocalNodeX3E1 - LocalNodeX2E1
g2E1 = LocalNodeX1E1 - LocalNodeX3E1
g3E1 = LocalNodeX2E1 - LocalNodeX1E1

StrainDisplacementMat1 = (1 / (2 * Area1)) * np.array([
    [ b1E1,      0, b2E1,      0,  b3E1,      0],
    [ 0,    g1E1,     0, g2E1,     0, g3E1],
    [g1E1, b2E1, g2E1,  b2E1, g3E1,  b3E1]
])

MaterialPropertyMatrix1 = (E / (1 - nu**2)) * np.array([
    [1,   nu,      0],
    [nu,   1,       0],
    [0,    0, ((1 - nu) / 2)]
])

k1 = MaterialThickness * Area1 * StrainDisplacementMat1.T @ MaterialPropertyMatrix1 @ StrainDisplacementMat1

# Element 2
LocalNodeX1E2 = 0.0
LocalNodeY1E2 = 0.0
LocalNodeX2E2 = RectLength
LocalNodeY2E2 = 0.0
LocalNodeX3E2 = RectLength
LocalNodeY3E2 = RectHeight

b1E2 = LocalNodeY2E2 - LocalNodeY3E2
b2E2 = LocalNodeY3E2 - LocalNodeY1E2
b3E2 = LocalNodeY1E2 - LocalNodeY2E2

# Note: Following the MATLAB code (which uses Element 1's g-values)
g1E2 = LocalNodeX3E2 - LocalNodeX2E2
g2E2 = LocalNodeX1E2 - LocalNodeX3E2
g3E2 = LocalNodeX2E2 - LocalNodeX1E2

StrainDisplacementMat2 = (1 / (2 * Area1)) * np.array([
    [b1E2,      0, b2E2,      0, b3E2,      0],
    [0,    g1E2,     0, g2E2,     0, g3E2],
    [g1E2, b2E2, g2E2, b2E2, g3E2, b3E2]
])

MaterialPropertyMatrix2 = (E / (1 - nu**2)) * np.array([
    [1,   nu,      0],
    [nu,  1,       0],
    [0,   0, (1 - nu) / 2]
])

k2 = MaterialThickness * Area1 * StrainDisplacementMat2.T @ MaterialPropertyMatrix2 @ StrainDisplacementMat2

# -------------------------------
# Global Stiffness Matrix Assembly
# -------------------------------
KGlobal = np.zeros((8, 8))


Elem1Dofs = [0, 1, 2, 3, 4, 5]
for i in range(6):
    for j in range(6):
        KGlobal[Elem1Dofs[i], Elem1Dofs[j]] += k1[i, j]


Elem2Dofs = [0, 1, 4, 5, 6, 7]
for i in range(6):
    for j in range(6):
        KGlobal[Elem2Dofs[i], Elem2Dofs[j]] += k2[i, j]

FixedDofs = [0, 1]  
FreeDofs = [i for i in range(8) if i not in FixedDofs]

KReduced = KGlobal[np.ix_(FreeDofs, FreeDofs)]
FGlobal = np.zeros((8, 1))
FGlobal[3, 0] = -Force 

DisplacementVect = np.zeros((8, 1))
DisplacementVect_free = np.linalg.solve(KReduced, FGlobal[FreeDofs])
for idx, dof in enumerate(FreeDofs):
    DisplacementVect[dof, 0] = DisplacementVect_free[idx]

# -------------------------------
# Compute Element Strains and Stresses for Element 1
# -------------------------------
Element1Dofs = [0, 1, 2, 3, 4, 5]
U1 = DisplacementVect[np.ix_(Element1Dofs)].flatten()

Epsilon1 = np.abs(StrainDisplacementMat1 @ U1)
Sigma1 = np.abs(MaterialPropertyMatrix1 @ Epsilon1)

print('Strains for element 1:')
print(Epsilon1)
print('Stresses for element 1:')
print(Sigma1)

# -------------------------------
# Compute Element Strains and Stresses for Element 2
# -------------------------------
Element2Dofs = [0, 1, 6, 7, 4, 5]
U2 = DisplacementVect[np.ix_(Element2Dofs)].flatten()

Epsilon2 = StrainDisplacementMat2 @ U2
Sigma2 = MaterialPropertyMatrix2 @ Epsilon2

print('Strains for element 2:')
print(Epsilon2)
print('Stresses for element 2:')
print(Sigma2)

# -------------------------------
# Stress-Strain Plot with Polynomial Fit for Both Elements
# -------------------------------
PolynomialDegree = 2

PolynomialFit1 = np.polyfit(Epsilon1, Sigma1, PolynomialDegree)
PolynomialFit2 = np.polyfit(Epsilon2, Sigma2, PolynomialDegree)

EpsilonFit1 = np.linspace(np.min(Epsilon1), np.max(Epsilon1), 100)
SigmaFit1 = np.polyval(PolynomialFit1, EpsilonFit1)

EpsilonFit2 = np.linspace(np.min(Epsilon2), np.max(Epsilon2), 100)
SigmaFit2 = np.polyval(PolynomialFit2, EpsilonFit2)

plt.figure()
plt.scatter(Epsilon1, Sigma1, color='red', label='Element 1 Data')
plt.scatter(Epsilon2, Sigma2, color='blue', label='Element 2 Data')
plt.plot(EpsilonFit1, SigmaFit1, '--', linewidth=2, label=f'{PolynomialDegree}th Order Fit')
plt.plot(EpsilonFit2, SigmaFit2, '--', linewidth=2, label=f'{PolynomialDegree}th Order Fit')
plt.xlabel('Strain')
plt.ylabel('Stress (Pa)')
plt.title(f'{PolynomialDegree}th Order Polynomial Fit for Stress-Strain Data')
plt.legend()
plt.grid(True)
plt.show()

print('Polynomial Fit Coefficients for Element 1:')
print(PolynomialFit1)

# -------------------------------
# Heatmap for U1 (X-displacement) and U2 (Y-displacement)
# -------------------------------
XNodes = np.array([0, RectLength, 0, RectLength])
YNodes = np.array([0, 0, RectHeight, RectHeight])

U1_nodes = np.array([DisplacementVect[0, 0], DisplacementVect[2, 0], DisplacementVect[4, 0], DisplacementVect[6, 0]])

U2_nodes = np.array([DisplacementVect[1, 0], DisplacementVect[3, 0], DisplacementVect[5, 0], DisplacementVect[7, 0]])

Xq, Yq = np.meshgrid(np.linspace(0, RectLength, 50), np.linspace(0, RectHeight, 50))
U1_interp = griddata((XNodes, YNodes), U1_nodes, (Xq, Yq), method='cubic')
U2_interp = griddata((XNodes, YNodes), U2_nodes, (Xq, Yq), method='cubic')

plt.figure()
plt.contourf(Xq, Yq, U1_interp, levels=20, cmap='viridis')
plt.colorbar()
plt.title('Displacement U1 (X-direction)')
plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.show()

plt.figure()
plt.contourf(Xq, Yq, U2_interp, levels=20, cmap='viridis')
plt.colorbar()
plt.title('Displacement U2 (Y-direction)')
plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.show()
