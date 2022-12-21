%% Analyze the results of RBCMotionInUnboundedFlows
%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon] = GridOnSphere(N);

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
ErrMemForce = zeros(Numframe,1);
ErrMemMomentForce = zeros(Numframe,1);

Volume = zeros(Numframe,1);
Area = zeros(Numframe,1);

dA = zeros(Numframe,1);
dB = zeros(Numframe,1);
dC = zeros(Numframe,1);
DeformationIndex = zeros(Numframe,1);

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

    [F, M, IntegralnormF, IntegralnormM] = integrateForce(axi, bxi, f);
    ErrMemForce(nframe) = (norm(F)/IntegralnormF)*100;
    ErrMemMomentForce(nframe) = (norm(M)/IntegralnormM)*100;

    [Volume(nframe), Area(nframe)] = getVolumeArea(axi, bxi);

    [dA(nframe), dB(nframe), dC(nframe)] = getDeformationIndex(axi, bxi);
    DeformationIndex(nframe) = ((dA(nframe)/dA(1)) - (dB(nframe)/dB(1))) / ...
                               ((dA(nframe)/dA(1)) + (dB(nframe)/dB(1)));
end
%% Dimensionalization
T_step = T_step/RefShearRate; % in seconds
Area = Area*(RefLength*10^(6))^2; % \mum^2
Volume = Volume*(RefLength*10^(6))^3; % \mum^3
ErrorVolume = max(abs(Volume(1)-max(Volume))/(Volume(1)), ...
                  abs(Volume(1)-min(Volume))/(Volume(1)))*100
ErrorArea = max(abs(Area(1)-max(Area))/(Area(1)), ...
               abs(Area(1)-min(Area))/(Area(1)))*100

dA = dA * (RefLength*10^(6)); % \mum
dB = dB * (RefLength*10^(6)); % \mum
dC = dC * (RefLength*10^(6)); % \mum