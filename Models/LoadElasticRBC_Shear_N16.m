%% Load the model specific parameters
%% 
name = 'ElasticRBC_Shear_N16';

%% Reference parameters 
ReferenceParameters

%% Parameters 
%% Fluid properties
%% Viscosity of exterior fluid
mu_out = 1.2*10^(-3); % Pa.s 
mu = mu_out/RefViscosity; 
%% Viscosity of interior fluid
mu_in = 6*10^(-3); % Pa.s
lam = mu_in/RefViscosity; % lam*mu

%% Shell properties
%% Material properties 
ShearModulus = 5.3*10^(-6); % N/m
ES = ShearModulus/RefElasticModulus; % inplane shear modulus
BendingModulus = 2*10^(-19); % N.m
EB = BendingModulus/RefBendingModulus; % bending modulus
DilatationModulus = 50*ShearModulus; % N/m
ED = DilatationModulus/RefElasticModulus; % inplane dilatational modulus

%% Standard-Linear-Solid SLS model parameters
%% Membrane viscosity
MembraneViscoelasticity = false;
mu_Membrane = 3.18*10^(-7); % m.Pa.s=N.s/m
mu_Mem = mu_Membrane/(RefLength*RefViscosity); % Bq, the Boussinesq number
ArtificialSpring = 100/3 * ShearModulus; % N/m 
% ArtificialSpring = 6*ShearModulus; % N/m 
RelaxationTime = mu_Membrane / ArtificialSpring; % s 
Tau = RelaxationTime*RefShearRate; % the relaxation time of the Maxwell element

%% Constitutive equations for membrane mechanics
StVenantKirchhoff = false;
NeoHookean = false;
Skalak = true;

%% Degree of spherical harmonic expansion
N = 16;
UpSampleFactor = 2; % Up-sample ratio for the force calculations
NGSphere = 2*N; % Quadrature on sphere parameter

%% Initial position and orientation
InitXi = [0; 0; 0];
InitOrient = expm(hat([0; 0; 0]*pi/2));

%% The tolerance of the GMRES method
ToleranceGMRES = 10^(-6);

%% Time step
%% Forward Euler
TimeIncrement = 10^(-5); % s
DT = TimeIncrement*RefShearRate;
EndTime_ = 0.15; % s
EndTime = EndTime_*RefShearRate;
NSTEPS = ceil(EndTime/DT);

SaveAtIncrementalSteps = 50; % Save the data at the increment of this steps

%% Ambient flow field
%% Simple shear flow
ShearFlow = true;
ShearRate_ = 123; % 1/s
ShearRate = ShearRate_/RefShearRate;
CapillaryNumberShear = ShearRate_ * mu_out * RefLength / ShearModulus;

%% Parabolic flow to mimic Poiseuille flow in a vessel
ParabolicFlow = false;
aParabolic_ = 1/(RefLength/RefShearRate); % 1/(m.s)
% For a vessel radius of 5*10^(-6) m 
bParabolic_ = (5*10^(-6))^2; % m^2
aParabolic = aParabolic_ * (RefLength/RefShearRate);
bParabolic = bParabolic_ / RefLength^2;
CapillaryNumberParabolic = mu_out * aParabolic_ * RefLength^4 / ...
                                                            BendingModulus;

%% Relaxation 
Relaxation = false;
RelaxationStep = NSTEPS; % relative to end time