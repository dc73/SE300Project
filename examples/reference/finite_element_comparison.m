clc
clear
close all

%% Inputs section

nu = 0.33;
E = 410; %GPa


RectHeight = input("Input the height of the rectangle structure [mm]: ");
RectLength = input("Input the length of the rectangle structure [mm]: ");
Force = input("Input single force's magnitude in [N]: ");
MaterialThickness = input("Input the thickness of the material [mm]: ");


ForceXComp = input("The X coordinate of the force on the beam: ");


while (ForceXComp < 0 || ForceXComp > RectLength)

    fprintf("Invalid value! Please enter a value between 0 and %.2f.\n", RectLength);

    ForceXComp = input("The X coordinate of the force on the beam: ");
end

fprintf("The X component of the force is logged.\n");

%% Calcualtion Section and Stiffness Matrix Calculation
% [k] = t A [B]^T [D] [B]

Area1 = 0.5*(RectHeight)*(RectLength);

%% Element 1
LocalNodeX1E1 = 0;
LocalNodeY1E1 = 0;
LocalNodeX2E1 = RectLength;
LocalNodeY2E1 = RectHeight;
LocalNodeX3E1 = 0;
LocalNodeY3E1 = RectHeight;

b1E1 = LocalNodeY2E1 - LocalNodeY3E1;
b2E1 = LocalNodeY3E1 - LocalNodeY1E1;
b3E1 = LocalNodeY1E1 - LocalNodeY2E1;

g1E1 = LocalNodeX3E1 - LocalNodeX2E1;
g2E1 = LocalNodeX1E1 - LocalNodeX3E1;
g3E1 = LocalNodeX2E1 - LocalNodeX1E1;

StrainDisplacementMat1= (1/(2*Area1))*[b1E1 0 b2E1 0 b3E1 0;
    0 g1E1 0 g2E1 0 g3E1;
    g1E1 b2E1 g2E1 b2E1 g3E1 b3E1];

MaterialPropertyMatrix1 = E/(1-nu^2)*[1 nu 0;
    nu 1 0;
    0 0 (1-nu)/2];

k1 = MaterialThickness*Area1*transpose(StrainDisplacementMat1)*MaterialPropertyMatrix1*StrainDisplacementMat1;

%% Element 2
LocalNodeX1E2 = 0;
LocalNodeY1E2 = 0;
LocalNodeX2E2 = RectLength;
LocalNodeY2E2 = 0;
LocalNodeX3E2 = RectLength;
LocalNodeY3E2 = RectHeight;

b1E2 = LocalNodeY2E2 - LocalNodeY3E2;
b2E2 = LocalNodeY3E2 - LocalNodeY1E2;
b3E2 = LocalNodeY1E2 - LocalNodeY2E2;

g1E2 = LocalNodeX3E1 - LocalNodeX2E1;
g2E2 = LocalNodeX1E1 - LocalNodeX3E1;
g3E2 = LocalNodeX2E1 - LocalNodeX1E1;

StrainDisplacementMat2= (1/(2*Area1))*[b1E2 0 b2E2 0 b3E2 0;
    0 g1E2 0 g2E2 0 g3E2;
    g1E2 b2E2 g2E2 b2E2 g3E2 b3E2];

MaterialPropertyMatrix2 = E/(1-nu^2)*[1 nu 0;
    nu 1 0;
    0 0 (1-nu)/2];

k2 = MaterialThickness*Area1*transpose(StrainDisplacementMat2)*MaterialPropertyMatrix2*StrainDisplacementMat2;

%% Defining the global stiffness matrix
KGlobal = zeros(8,8);

Elem1Dofs = [1,2, 3,4, 5,6];
for i = 1:6
    for j = 1:6
        KGlobal(Elem1Dofs(i), Elem1Dofs(j)) = ...
            KGlobal(Elem1Dofs(i), Elem1Dofs(j)) + k1(i,j);
    end
end

Elem2Dofs = [1,2, 5,6, 7,8];
for i = 1:6
    for j = 1:6
        KGlobal(Elem2Dofs(i), Elem2Dofs(j)) = ...
            KGlobal(Elem2Dofs(i), Elem2Dofs(j)) + k2(i,j);
    end
end

FixedDofs = [1,2];
FreeDofs = setdiff(1:8, FixedDofs);

KReduced = KGlobal(FreeDofs, FreeDofs);
FGlobal = zeros(8,1);
FGlobal(4) = -Force;

DisplacementVect = zeros(8,1);
DisplacementVect(FreeDofs) = KReduced \ FGlobal(FreeDofs);

%% Compute Element Strains and Stresses for Element 1
Element1Dofs = [1, 2, 3, 4, 5, 6];
U1 = DisplacementVect(Element1Dofs);

Epsilon1 = abs(StrainDisplacementMat1 * U1);
Sigma1 = abs(MaterialPropertyMatrix1 * Epsilon1);

disp('Strains for element 1:');
disp(Epsilon1);
disp('Stresses for element 1:');
disp(Sigma1);

%% Compute Element Strains and Stresses for Element 2
Element2Dofs = [1, 2, 7, 8, 5, 6];
U2 = DisplacementVect(Element2Dofs);

Epsilon2 = StrainDisplacementMat2 * U2;
Sigma2 = MaterialPropertyMatrix2 * Epsilon2;

disp('Strains for element 2:');
disp(Epsilon2);
disp('Stresses for element 2:');
disp(Sigma2);

%% Stress-Strain Plot
PolynomialDegree = 3; % Change this value to 2, 3, or 4 as needed
PolynomialFit1 = polyfit(Epsilon1, Sigma1, PolynomialDegree);
PolynomialFit2 = polyfit(Epsilon2, Sigma2, PolynomialDegree);

% Generate smooth curve data for plotting
EpsilonFit1 = linspace(min(Epsilon1), max(Epsilon1), 100);
SigmaFit1 = polyval(PolynomialFit1, EpsilonFit1);

EpsilonFit2 = linspace(min(Epsilon2), max(Epsilon2), 100);
SigmaFit2 = polyval(PolynomialFit2, EpsilonFit2);

% Plot original data and fitted curve
figure;
scatter(Epsilon1, Sigma1, 'ro', 'filled'); % Original data points
scatter(Epsilon2,Sigma2,'bo', 'filled')
hold on;
plot(EpsilonFit1, SigmaFit1, '--', 'LineWidth', 2, 'DisplayName', sprintf('%dth Order Fit', PolynomialDegree));
plot(EpsilonFit2, SigmaFit2, '--', 'LineWidth', 2, 'DisplayName', sprintf('%dth Order Fit', PolynomialDegree));
hold off;
hold off;

xlabel('Strain');
ylabel('Stress (Pa)');
title(sprintf('%dth Order Polynomial Fit for Stress-Strain Data', PolynomialDegree));
legend;
grid on;

figure;
scatter(Epsilon1, Sigma1, 'ro', 'filled'); % Element 1 data points
hold on;
scatter(Epsilon2, Sigma2, 'bo', 'filled'); % Element 2 data points
plot(EpsilonFit1, SigmaFit1, '--', 'LineWidth', 2, 'DisplayName', 'Element 1 Fit');
plot(EpsilonFit2, SigmaFit2, '-.', 'LineWidth', 2, 'DisplayName', 'Element 2 Fit');
hold off;

xlabel('Strain');
ylabel('Stress (Pa)');
title(sprintf('%dth Order Polynomial Fit for Combined Stress-Strain Data', PolynomialDegree));
legend;
grid on;


% Display polynomial coefficients
disp('Polynomial Fit Coefficients:');
disp(PolynomialFit1);



%% Heatmap for U1 (X-displacement) and U2 (Y-displacement)
XNodes = [0, RectLength, 0, RectLength];
YNodes = [0, 0, RectHeight, RectHeight];

U1 = [DisplacementVect(1), DisplacementVect(3), DisplacementVect(5), DisplacementVect(7)];
U2 = [DisplacementVect(2), DisplacementVect(4), DisplacementVect(6), DisplacementVect(8)];

[Xq, Yq] = meshgrid(linspace(0, RectLength, 50), linspace(0, RectHeight, 50));

U1_interp = griddata(XNodes, YNodes, U1, Xq, Yq, 'cubic');
U2_interp = griddata(XNodes, YNodes, U2, Xq, Yq, 'cubic');

figure;
contourf(Xq, Yq, U1_interp, 20, 'LineColor', 'none');
colorbar;
title('Displacement U1 (X-direction)');
xlabel('X [mm]');
ylabel('Y [mm]');

figure;
contourf(Xq, Yq, U2_interp, 20, 'LineColor', 'none');
colorbar;
title('Displacement U2 (Y-direction)');
xlabel('X [mm]');
ylabel('Y [mm]');
