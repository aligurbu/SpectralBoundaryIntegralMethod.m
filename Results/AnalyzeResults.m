%% Analyze the results
addpath(genpath('../../SpectralBoundaryIntegralMethod.m'))

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16
% LoadElasticRBC_Parabolic_N16
% LoadMemViscosityRBC_Shear_N16
% LoadMemViscosityRBC_Parabolic_N16

%% Set up output file
fidCoord = fopen(['Coord_',name,'.dat'],'r');
fidMemFor = fopen(['MemFor_',name,'.dat'],'r');
fidSol = fopen(['Sol_',name,'.dat'],'r');
Time = fread(fidTime,'double');
NSTEPS = length(Time);

