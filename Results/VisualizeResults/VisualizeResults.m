%% Visualize the results
clear all; close all; clc;
addpath(genpath('../../../SpectralBoundaryIntegralMethod.m'))
verbose_Plot = false;
WritetoGIF = true; % if false then save as MPEG-4

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16
% LoadElasticRBC_Parabolic_N16
% LoadMemViscosityRBC_Shear_N16
% LoadMemViscosityRBC_Parabolic_N16

%% Set up output file
fidCoord = fopen(['Coord_',name,'.dat'],'r');
fidMemFor = fopen(['MemFor_',name,'.dat'],'r');
fidSol = fopen(['Sol_',name,'.dat'],'r');
fidTime = fopen(['Time_',name,'.dat'],'r');
Time = fread(fidTime,'double');
NSTEPS = length(Time);

%%
timeStepIncrement = 1;
FRAMES = floor(linspace(1,NSTEPS,5));
TimesOfFrames = Time(FRAMES)/RefShearRate % in seconds

Rendering = 2; % This upsampling for the rendering purposes

%% Choose a view angle
% viewInd = [-45 30];
% name = [name,'_3D'];
viewInd = [0 90];
name = [name,'_xy'];

%% Choose a background color for visualization
blackBackground = true; % if false then white background

%% Run Visualize m-files one at a time.
VisualizeMembraneShape

% CaptureFramesMembraneShape

% VisualizeIsotropicMembraneTension

% CaptureFramesIsotropicMembraneTension

% VisualizeMembraneForceVectorField

% VisualizeMembraneVelocityVectorField

fclose('all');