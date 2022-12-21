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
