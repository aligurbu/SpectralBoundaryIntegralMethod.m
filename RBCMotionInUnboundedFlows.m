%% Spectral Galerkin BIE implementation for an RBC motion and deformation
%% in an unbounded domain under the ambient flow field
%%
clear all; close all; clc;
addpath(genpath('../SpectralBoundaryIntegralMethod.m'))

verbose_Plot = false;
Starttime = tic;

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

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%% Set up output file
fidTime = fopen(['Time_',name,'.dat'],'w');
fidCoord = fopen(['Coord_',name,'.dat'],'w');
fidMemFor = fopen(['MemFor_',name,'.dat'],'w');
fidSol = fopen(['Sol_',name,'.dat'],'w');

%% Set up time integration
cu = zeros(3*(N+1)^2,1);

%% Initialization
xi = Xi; axi = aXi; bxi = bXi;

viscousStress_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);
epsilbrev_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);

%% Time-stepping
for nstep = 0:1%NSTEPS
    if (nstep==0 || mod(nstep,SaveAtIncrementalSteps)==0)
        %% Write the time to file
        fwrite(fidTime, Time, 'double');
        %% Write the position of cell to file
        cxi = [axi(mask_a); bxi(mask_b)];
        fwrite(fidCoord, cxi, 'double');
    end

    %% Ambient simple shear flow
    if ShearFlow
        uinf = zeros(size(Xi));
        uinf(:,:,1) = ShearRate*xi(:,:,2); % Shear in y-direction
%         uinf(:,:,1) = ShearRate*xi(:,:,3); % Shear in z-direction
        [auinf,buinf] = shagcm(uinf);
        if Relaxation && nstep > RelaxationStep
            uinf = zeros(size(Xi));
            auinf = zeros(size(axi));
            buinf = zeros(size(bxi));
        end
    end
    %% Ambient parabolic velocity field
    if ParabolicFlow
        uinf = zeros(size(Xi));
        uinf(:,:,1) = aParabolic*(bParabolic-(xi(:,:,2).^2 + xi(:,:,3).^2));
        [auinf,buinf] = shagcm(uinf);
        if Relaxation && nstep > RelaxationStep
            uinf = zeros(size(Xi));
            auinf = zeros(size(axi));
            buinf = zeros(size(bxi));
        end
    end


    %% Compute membrane traction
    [f, viscousStress_prev, epsilbrev_prev] = ...
                           MembraneForces(EE, FF, GG, WW, JXibrev, ...
                                          LL, MM, NN, ...
                                          axi, bxi, ...
                                          UpSampleFactor, ...
                                          ES, ED, EB, ...
                                          StVenantKirchhoff, ...
                                          NeoHookean, Skalak, ...
                                          MembraneViscoelasticity, ...
                                          mu_Mem, Tau, DT, ...
                                          viscousStress_prev, ...
                                          epsilbrev_prev);
    [af,bf] = shagcm(f);

    %% Solve Boundary Integral Equation
    [u, au, bu, cu, FLAG, RELRES, ITER] = ...
                           SolveSpectralBIE(axi, bxi, xi, ...
                                            RotationMatrix, ...
                                            NGSphere, NGthet, NGphi, ...
                                            eta, wg, ...
                                            af, bf, auinf, buinf, ...
                                            nlat, nlon, mask_a, mask_b, ...
                                            mu, lam, N, ToleranceGMRES, ...
                                            cu);
    if (nstep==0 || mod(nstep,SaveAtIncrementalSteps)==0)
        %% Write output to file
        writef = f(:)';
        fwrite(fidMemFor,writef,'double');
        writecu = cu';
        fwrite(fidSol,writecu,'double');
    end

	%% Update state (Forward Euler time scheme)
    axi = axi + DT*au;
    bxi = bxi + DT*bu;
    xi = shsgcm(axi,bxi);

    %% Set the center of the mass of RBC to the initial position, InitXi
    if ParabolicFlow
        [gxi_thet, gxi_phi] = gradgcm(axi, bxi);
        normalVector = cross(gxi_thet, gxi_phi, 3);
        xi_dot_normalVector = sum(xi.*normalVector,3);
        xi_xi_dot_normalVector = xi.*xi_dot_normalVector;
        V = sum(sum(xi_dot_normalVector.*weight))*(2*pi/nlon)/3;
        Int_xi_xi_dot_normalVector_dA = ...
                    sum(sum(xi_xi_dot_normalVector.*weight))*(2*pi/nlon)/4;
        CenterMass = reshape(Int_xi_xi_dot_normalVector_dA./V,3,1);
        xi = xi - reshape((CenterMass - InitXi),1,1,3);
    end
    [axi,bxi] = shagcm(xi);






































end
fclose('all');
RunTimeInHours = toc(Starttime)/60/60