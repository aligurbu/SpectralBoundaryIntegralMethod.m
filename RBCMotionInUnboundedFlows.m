%% Spectral Galerkin BIE implementation for an RBC motion and deformation
%% in an unbounded domain under the ambient flow field
%%
clear all; close all; clc;
addpath(genpath('../SpectralBoundaryIntegralMethod.m'))

verbose_Plot = false;

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16

%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon, thet, phi, weight] = GridOnSphere(N);

%% Create grid for integration points on unit sphere:
[NGthet, NGphi, eta, wg] = setup_integration_grid(N, NGSphere);

%% Create geometry of RBC: undeformed geometry of an RBC
Xi = getRBCInitialGeometry(thet,phi,InitXi,InitOrient);

%% Initial conditions:
%% spherical harmonics coefficients of the undeformed RBC coordinates
[aXi,bXi] = shagcm(Xi);

%% Compute the first and second fundamental form coefficients
[~, ~, EE, FF, GG, WW, JXibrev, ~, ~, ~, LL, MM, NN] = ...
                   coefficientsOfFundamentalForm(aXi, bXi, UpSampleFactor);

%% Check the geometry and position of RBC
if verbose_Plot
    figure('Color','white')
    hold on
    VisualizeGeometry(nlat, nlon, aXi, bXi, 'r', true)
    axis on
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

%% Create rotation matrices
RotationMatrix = RotationMatrixForSHCoefficients(N, nlat, nlon, thet, phi);