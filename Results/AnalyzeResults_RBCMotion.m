%% Analyze the results of RBCMotionInUnboundedFlows
%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon, thet, phi, weight] = GridOnSphere(N);

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%%
Numframe = 0;
for nstep = 0:NSTEPS-1
    if (nstep~=0 && nstep~=NSTEPS-1 && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    Numframe = Numframe + 1;
end
%%
viscousStress_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);
epsilbrev_prev = zeros(UpSampleFactor*N+1, 2*UpSampleFactor*N+1, 4);

%%
F_ = zeros(3,Numframe);
ErrMemForce_ = zeros(Numframe,1);
M_ = zeros(3,Numframe);
ErrMemMomentForce_ = zeros(Numframe,1);

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 0:NSTEPS-1
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    f = reshape(fread(fidMemFor,nlat*nlon*3,'double'),nlat,nlon,3);

    aU = fread(fidSol,3*(N+1)*(N+2)/2,'double');
    bU = fread(fidSol,3*N*(N+1)/2,'double');

    if (nstep~=0 && nstep~=NSTEPS-1 && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep+1);

    au = zeros(size(mask_a)); bu = zeros(size(mask_b));
    au(mask_a) = aU; bu(mask_b) = bU;
    u = shsgcm(au, bu);

    %%
    [F__, M__, IntnormF, IntnormM] = integrateForce(axi, bxi, f);
    F_(:,nframe) = F__;
    M_(:,nframe) = M__;
    ErrMemForce_(nframe) = (norm(F__)/IntnormF)*100;
    ErrMemMomentForce_(nframe) = (norm(M__)/IntnormM)*100;
end
%% Dimensionalization
T_step = T_step/RefShearRate; % in seconds