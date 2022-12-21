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