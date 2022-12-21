%% Spectral Galerkin BIE implementation for an RBC motion and deformation
%% in an unbounded domain under the ambient flow field
%%
clear all; close all; clc;
addpath(genpath('../SpectralBoundaryIntegralMethod.m'))

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16

%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon, thet, phi, weight] = GridOnSphere(N);

%% Create grid for integration points on unit sphere:
[NGthet, NGphi, eta, wg] = setup_integration_grid(N, NGSphere);

%% Create geometry of RBC: undeformed geometry of an RBC
Xi = getRBCInitialGeometry(thet,phi,InitXi,InitOrient);