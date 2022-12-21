%% Post-processing the flow around the RBC motion

%% Create grid for representation: [N+1 x 2*N+1]
[nlat, nlon, thet, phi, ~] = GridOnSphere(N);

%% Create fine equally-spaced grid for nearly-singular integrals
[nlatEqSpaced, nlonEqSpaced, phiEqSpaced, thetEqSpaced] = ...
                                     EquallySpacedGridOnSphere(N_EqSpaced);

%% Create rotation matrices
RotationMatrix_EqSpaced = ...
         RotationMatrixForSHCoefficients(N, nlatEqSpaced, nlonEqSpaced, ...
                                         phiEqSpaced, thetEqSpaced);

%% Create grid for integration points on unit sphere:
%% Note that eta is NOT needed for the nearly-singular integration!
[NGthet, NGphi, ~, wg] = setup_integration_grid(N, NGSphere);

%% Create geometry of RBC: undeformed geometry of an RBC
Xi = getRBCInitialGeometry(thet,phi,InitXi,InitOrient);
% spherical harmonics coefficients of the undeformed RBC coordinates
[aXi,bXi] = shagcm(Xi);

% Check the geometry and position of RBC and depict EvaluationPts
if verbose_Plot
    figure('Color','white')
    hold on
    VisualizeGeometry(nlat, nlon, aXi, bXi, 'r', true)
    plot3(EvaluationPts(1,:), EvaluationPts(2,:), EvaluationPts(3,:), ...
          'k.','MarkerSize',MarkerSizeind)
    axis on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(viewInd)
end

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

%%
Numframe = 0;
for nstep = 1:NSTEPS
    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    Numframe = Numframe + 1;
end

%%
xi_GaussX = zeros(nlat,nlon,Numframe);
xi_GaussY = zeros(nlat,nlon,Numframe);
xi_GaussZ = zeros(nlat,nlon,Numframe);
CellX = zeros(size(faces,1),size(faces,2),Numframe);
CellY = zeros(size(faces,1),size(faces,2),Numframe);
CellZ = zeros(size(faces,1),size(faces,2),Numframe);

u_GaussX = zeros(nlat,nlon,Numframe);
u_GaussY = zeros(nlat,nlon,Numframe);
u_GaussZ = zeros(nlat,nlon,Numframe);

VelEvalPtsX = zeros(numPoints,Numframe);
VelEvalPtsY = zeros(numPoints,Numframe);
VelEvalPtsZ = zeros(numPoints,Numframe);
SpeedEvalPts = zeros(numPoints,Numframe);

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 1:NSTEPS
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    f = reshape(fread(fidMemFor,nlat*nlon*3,'double'),nlat,nlon,3);

    aU = fread(fidSol,3*(N+1)*(N+2)/2,'double');
    bU = fread(fidSol,3*N*(N+1)/2,'double');

    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep);

    %% Finer equally-spaced grid for the closest point computation
    axiEqSpaced = UpSampling(axi,N,N_EqSpaced);
    bxiEqSpaced = UpSampling(bxi,N,N_EqSpaced);
    xiEqSpaced_ = shsecm(axiEqSpaced,bxiEqSpaced);
    xiEqSpaced = reshape(xiEqSpaced_,nlatEqSpaced*nlonEqSpaced,3)';

    %% Determine whether the target point, chi, is inside or outside
    % All the points outside for now.
    EvaluationPtsLocation = false(size(EvaluationPts,2),1);
    xi_Gauss = shsgcm(axi,bxi);
    xiGauss = reshape(xi_Gauss,nlat*nlon,3)';
    [gxi_thet,gxi_phi] = gradgcm(axi,bxi);
    Normal(:,:,1) = gxi_thet(:,:,2).*gxi_phi(:,:,3) - ...
                    gxi_thet(:,:,3).*gxi_phi(:,:,2);
    Normal(:,:,2) = gxi_thet(:,:,3).*gxi_phi(:,:,1) - ...
                    gxi_thet(:,:,1).*gxi_phi(:,:,3);
    Normal(:,:,3) = gxi_thet(:,:,1).*gxi_phi(:,:,2) - ...
                    gxi_thet(:,:,2).*gxi_phi(:,:,1);
    NormalGauss = reshape(Normal,nlat*nlon,3)';
    for n = 1 : size(EvaluationPts,2)
        chi = EvaluationPts(:,n);
        [~,IndGaussmin] = min(sqrt(sum((xiGauss-chi).^2,1)));
        XiGaussPts = xiGauss(:,IndGaussmin);

        NormalGaussPts = NormalGauss(:,IndGaussmin);
        if (chi - XiGaussPts)'*NormalGaussPts < 0
            % The target point inside
            EvaluationPtsLocation(n) = true;
        end
    end

    %% Ambient parabolic velocity field
    if ParabolicFlow
        uAMB = zeros(size(EvaluationPts));
        uAMB(1,:) = aParabolic*...
              (bParabolic-(EvaluationPts(2,:).^2 + EvaluationPts(3,:).^2));
        if Relaxation && nstep > RelaxationStep
            uAMB = zeros(size(EvaluationPts));
        end
    end
    %% Ambient simple shear flow
    if ShearFlow
        uAMB = zeros(size(EvaluationPts));
        uAMB(1,:) = ShearRate*EvaluationPts(2,:); % Shear in y-direction
%         uAMB(1,:) = ShearRate*EvaluationPts(3,:); % Shear in z-direction
        if Relaxation && nstep > RelaxationStep
            uAMB = zeros(size(EvaluationPts));
        end
    end
    uAMB(:,EvaluationPtsLocation) = uAMB(:,EvaluationPtsLocation)/lam;

    %% Compute the nearly-singular integral with Gf
    [af,bf] = shagcm(f);
    cf = [af(mask_a); bf(mask_b)];
    cxi = [axi(mask_a); bxi(mask_b)];
    Gf = NearlySingular_ComputeGf(cf, cxi, xiEqSpaced, ...
                                  EvaluationPts, EvaluationPtsLocation, ...
                                  RotationMatrix_EqSpaced, ...
                                  wg, mask_a, mask_b, mu, lam, ...
                                  N, NGSphere, NGphi);

	%% Compute the integrand for the nearly-singular integral with K kernel
    Kernel = NearlySingular_ComputeK(cxi, xiEqSpaced, ...
                                     EvaluationPts, ...
                                     RotationMatrix_EqSpaced, ...
                                     wg, mask_a, mask_b, ...
                                     N, NGSphere, NGthet, NGphi);

	%% Compute the nearly-singular integral with K(u-u_\perp)
    au = zeros(size(axi)); bu = zeros(size(bxi));
    au(mask_a) = aU; bu(mask_b) = bU;
    cu = [aU;bU];

    auEqSpaced = UpSampling(au,N,N_EqSpaced);
    buEqSpaced = UpSampling(bu,N,N_EqSpaced);
    uEqSpaced_ = shsecm(auEqSpaced,buEqSpaced);
    uEqSpaced = reshape(uEqSpaced_,nlatEqSpaced*nlonEqSpaced,3)';

    Ku = NearlySingular_ComputeKu(Kernel, cu, xiEqSpaced, uEqSpaced, ...
                                  EvaluationPts, EvaluationPtsLocation, ...
                                  RotationMatrix_EqSpaced, ...
                                  mask_a, mask_b, lam, ...
                                  N, NGSphere);

	%% Velocity at the evaluation points
    VelocityEvaluationPts = uAMB(:) - Ku - Gf;

    VelEvalPts = reshape(VelocityEvaluationPts,numDofPerNode,numPoints);

    VelEvalPtsX(:,nframe) = VelEvalPts(1,:);
    VelEvalPtsY(:,nframe) = VelEvalPts(2,:);
    VelEvalPtsZ(:,nframe) = VelEvalPts(3,:);

    SpeedEvalPts(:,nframe) = sqrt(VelEvalPts(1,:).^2 + ...
                                  VelEvalPts(2,:).^2 + ...
                                  VelEvalPts(3,:).^2);

    %% Coordinates of vertices
    xi_equal = shsecm(axi,bxi);
    xi_Gauss = shsgcm(axi,bxi);
    xi_GaussX(:,:,nframe) = xi_Gauss(:,:,1);
    xi_GaussY(:,:,nframe) = xi_Gauss(:,:,2);
    xi_GaussZ(:,:,nframe) = xi_Gauss(:,:,3);
    xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
    Vert = reshape(xiplot,(nlat+2)*nlon,3);
    CellX(:,:,nframe) = reshape(Vert(faces(:),1),size(faces));
    CellY(:,:,nframe) = reshape(Vert(faces(:),2),size(faces));
    CellZ(:,:,nframe) = reshape(Vert(faces(:),3),size(faces));

    %% Membrane velocity on patches vertices
    u_Gauss = shsgcm(au, bu);
    u_GaussX(:,:,nframe) = u_Gauss(:,:,1);
    u_GaussY(:,:,nframe) = u_Gauss(:,:,2);
    u_GaussZ(:,:,nframe) = u_Gauss(:,:,3);
end
%% Dimensionalization
T_step = T_step/RefShearRate; % in seconds
EvaluationPts = EvaluationPts*RefLength*10^(6); % \mum
EvaluationPtsX = EvaluationPtsX*RefLength*10^(6); % \mum
EvaluationPtsY = EvaluationPtsY*RefLength*10^(6); % \mum
CellX = CellX*RefLength*10^(6); % \mum
CellY = CellY*RefLength*10^(6); % \mum
CellZ = CellZ*RefLength*10^(6); % \mum
SpeedEvalPts = SpeedEvalPts*RefVelocity*10^(3); % mm/s
VelEvalPtsX = VelEvalPtsX*RefVelocity*10^(3); % mm/s
VelEvalPtsY = VelEvalPtsY*RefVelocity*10^(3); % mm/s
VelEvalPtsZ = VelEvalPtsZ*RefVelocity*10^(3); % mm/s

%%
if blackBackground
    figure('Color','black');
    ColorInd = 'white';
else
    figure('Color','white');
    ColorInd = 'black';
end
hold on
scr_siz = get(0,'ScreenSize');
scr_siz([1 2]) = scr_siz([1 2]) + 100;
scr_siz([3 4]) = scr_siz([3 4]) - 200;
set(gcf, 'Position',  scr_siz)

p1 = plot3(CellX(1,1,1), ...
           CellY(1,1,1), ...
           CellZ(1,1,1), ...
           'k.','MarkerSize',5*MarkerSizeind);
p2 = plot3(CellX(1,9,1), ...
           CellY(1,9,1), ...
           CellZ(1,9,1), ...
           'r.','MarkerSize',5*MarkerSizeind);
p3 = plot3(CellX(1,298,1), ...
           CellY(1,298,1), ...
           CellZ(1,298,1), ...
           'c.','MarkerSize',5*MarkerSizeind);

%
uSpeed = reshape(SpeedEvalPts(:,1),size(EvaluationPtsX));
% ax = gca;
% HG = hgtransform(ax);
[j, c] = contourf(EvaluationPtsX, EvaluationPtsY, uSpeed, 200, ...
                                        'LineColor', 'none');%, 'Parent', HG);
% HG.Matrix = makehgtform('xrotate', pi/2);
caxis([min(min(SpeedEvalPts)) max(max(SpeedEvalPts))])
colormap(jet)
cbar = colorbar('southoutside');
set(cbar,'FontSize',12,'Color',ColorInd)
set(get(cbar,'title'),'string','Speed (mm/s)','Color',ColorInd);
set(gca,'FontName','cambria math','FontSize',12)

t = quiver3(EvaluationPts(1,:)', EvaluationPts(2,:)', EvaluationPts(3,:)', ...
            VelEvalPtsX(:,1), VelEvalPtsY(:,1), VelEvalPtsZ(:,1), 1.5,...
            ColorInd);
h = patch(CellX(:,:,1),CellY(:,:,1),CellZ(:,:,1),'r');
alpha(h, TransparencyInd) % to set transparency
set(gca,'FontName','cambria math','FontSize',12)
material Dull
axis off
view(viewInd)
set(gca,'DataAspectRatio',[1 1 1])
ht = title({sprintf('Time = %4.4f sec',T_step(1))},'Color',ColorInd);
set(ht,'FontName','cambria math','FontSize',12)
axis equal
axis off
box off

%
v = VideoWriter('PostProcessing.mp4','MPEG-4');
v.FrameRate = 10;
open(v);

% Loop through by changing XData and YData
for id = 1:length(T_step)
    %% Update graphics data. This is more efficient than recreating plots.
    uSpeed = reshape(SpeedEvalPts(:,id),size(EvaluationPtsX));
    set(c, 'ZData', uSpeed);
    caxis(gca, [min(min(uSpeed)) max(max(uSpeed))]);
    set(t, 'UData', VelEvalPtsX(:,id), ...
           'VData', VelEvalPtsY(:,id), ...
           'WData', VelEvalPtsZ(:,id), 'Color', ColorInd);
    set(h, 'XData', CellX(:,:,id), ...
           'YData', CellY(:,:,id), ...
           'ZData', CellZ(:,:,id), ...
           'FaceAlpha', TransparencyInd);
    set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
    
    set(p1, 'XData', CellX(1,1,id), ...
            'YData', CellY(1,1,id), ...
            'ZData', CellZ(1,1,id));
    set(p2, 'XData', CellX(1,9,id), ...
            'YData', CellY(1,9,id), ...
            'ZData', CellZ(1,9,id));
    set(p3, 'XData', CellX(1,298,id), ...
            'YData', CellY(1,298,id), ...
            'ZData', CellZ(1,298,id));

    %% Get frame as an image
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);