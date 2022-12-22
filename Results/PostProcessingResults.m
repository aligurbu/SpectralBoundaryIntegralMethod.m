%% Post-processing the flow around the RBC motion
clear all; close all; clc;
addpath(genpath('../../SpectralBoundaryIntegralMethod.m'))
verbose_Plot = false;
Starttime = tic;

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16
% LoadElasticRBC_Parabolic_N16
% LoadMemViscosityRBC_Shear_N16
% LoadMemViscosityRBC_Parabolic_N16

N_EqSpaced = 2*N;

%% Set up output file
fidCoord = fopen(['Coord_',name,'.dat'],'r');
fidMemFor = fopen(['MemFor_',name,'.dat'],'r');
fidSol = fopen(['Sol_',name,'.dat'],'r');
fidTime = fopen(['Time_',name,'.dat'],'r');
Time = fread(fidTime,'double');
NSTEPS = length(Time);

%%
timeStepIncrement = 2;

%% Choose a view angle
viewInd = [0 90];

%% Choose a background color for visualization
blackBackground = true; % if false then white background

%% Visualization settings
VisualizeSettings
TransparencyInd = 0;

%% Set-up the evaluation points
evaluationPoitnsGridSize = sqrt(bParabolic);
Xpoints = linspace(-1.5*evaluationPoitnsGridSize,1.5*evaluationPoitnsGridSize,45);
Ypoints = linspace(-evaluationPoitnsGridSize,evaluationPoitnsGridSize,30);
[EvaluationPtsX, EvaluationPtsY] = meshgrid(Xpoints, Ypoints);
EvaluationPts = zeros(3, numel(EvaluationPtsX));
EvaluationPts(1,:) = EvaluationPtsX(:)';
EvaluationPts(2,:) = EvaluationPtsY(:)';
numDofPerNode = size(EvaluationPts,1);
numPoints = size(EvaluationPts,2);

%% Post-processing the flow around the RBC motion
PostProcessing_RBCMotion

%%
fclose('all');
RunTime = toc(Starttime)