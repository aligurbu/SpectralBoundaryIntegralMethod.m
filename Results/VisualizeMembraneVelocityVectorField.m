%% Visualize velocity vector fields on the cell membrane
%% Create grid for representation: [N+1 x 2*N+1]
nlat = N+1;
nlon = 2*N+1;

%% Visualization settings
VisualizeSettings
TransparencyInd = 0.25*TransparencyInd;

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

%% Masks to go between Spherepack and vector representations of SH coeff
mask_a = repmat(triu(true(N+1),0),1,1,3);
mask_b = mask_a;
mask_b(1,:,:) = false;

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

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 1:NSTEPS
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    aU = fread(fidSol,3*(N+1)*(N+2)/2,'double');
    bU = fread(fidSol,3*N*(N+1)/2,'double');

    if (nstep~=1 && nstep~=NSTEPS && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep)/RefShearRate; % in seconds

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
    au = zeros(size(axi)); bu = zeros(size(bxi));
    au(mask_a) = aU; bu(mask_b) = bU;
    u_Gauss = shsgcm(au, bu);
    u_GaussX(:,:,nframe) = u_Gauss(:,:,1);
    u_GaussY(:,:,nframe) = u_Gauss(:,:,2);
    u_GaussZ(:,:,nframe) = u_Gauss(:,:,3);
end
MinMaxBounds = [min(min(min(CellX))) max(max(max(CellX))) ...
                min(min(min(CellY))) max(max(max(CellY))) ...
                min(min(min(CellZ))) max(max(max(CellZ)))];
MinMaxBounds = MinMaxBounds + [-1 1 -1 1 -1 1]*0.25;

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

p = scatter3(reshape(xi_GaussX(:,:,1),nlat*nlon,1), ...
             reshape(xi_GaussY(:,:,1),nlat*nlon,1), ...
             reshape(xi_GaussZ(:,:,1),nlat*nlon,1), [ColorInd,'.']);
q = quiver3(xi_GaussX(:,:,1), xi_GaussY(:,:,1), xi_GaussZ(:,:,1), ...
            u_GaussX(:,:,1), u_GaussY(:,:,1), u_GaussZ(:,:,1), ...
            'filled', ColorInd);
h = patch(CellX(:,:,1),CellY(:,:,1),CellZ(:,:,1),'r');
alpha(h, TransparencyInd) % to set transparency
% set(h, 'linestyle', 'none') % no lines showning element edges
set(gca,'FontName','cambria math','FontSize',12)
material Dull
axis off
view(viewInd)
set(gca,'DataAspectRatio',[1 1 1])
ht = title({sprintf('Time = %4.4f sec',T_step(1))},'Color',ColorInd);
set(ht,'FontName','cambria math','FontSize',12)
axis(MinMaxBounds)
axis equal
axis off
box off

%% Get figure size
if WritetoGIF
    frame = getframe(gcf);
    height = size(frame.cdata, 1);
    width = size(frame.cdata, 2);

    %% Preallocate data (for storing frame data)
    mov = zeros(height, width, 1, length(T_step), 'uint8');
    
    %% Loop through by changing XData and YData
    for id = 1:length(T_step)
        %% Update graphics data. This is more efficient than recreating plots.
        set(p, 'XData', reshape(xi_GaussX(:,:,id),nlat*nlon,1), ...
               'YData', reshape(xi_GaussY(:,:,id),nlat*nlon,1), ...
               'ZData', reshape(xi_GaussZ(:,:,id),nlat*nlon,1));
        set(q, 'XData', xi_GaussX(:,:,id), ...
               'YData', xi_GaussY(:,:,id), ...
               'ZData', xi_GaussZ(:,:,id), ...
               'UData', u_GaussX(:,:,id), ...
               'VData', u_GaussY(:,:,id), ...
               'WData', u_GaussZ(:,:,id), 'Color', ColorInd);
        set(h, 'XData', CellX(:,:,id), ...
               'YData', CellY(:,:,id), ...
               'ZData', CellZ(:,:,id), ...
               'FaceAlpha', TransparencyInd);%, 'linestyle', 'none')
        set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
        axis(MinMaxBounds)
        %% Get frame as an image
        frame = getframe(gcf);
    
        %% Create a colormap for the first frame. For the rest of the frames,
        %% use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(frame.cdata, 1024, 'nodither');
        else
            mov(:,:,1,id) = rgb2ind(frame.cdata, map, 'nodither');
        end
    end
    %% Create animated GIF
    imwrite(mov, map, 'MembraneVelocityProfile.gif', 'DelayTime', 0, 'LoopCount', inf)
else % Write to MP4 file
    v = VideoWriter('MembraneVelocityProfile.mp4','MPEG-4');
    v.FrameRate = 10;
    open(v);

    %% Loop through by changing XData and YData
    for id = 1:length(T_step)
        %% Update graphics data. This is more efficient than recreating plots.
        set(p, 'XData', reshape(xi_GaussX(:,:,id),nlat*nlon,1), ...
               'YData', reshape(xi_GaussY(:,:,id),nlat*nlon,1), ...
               'ZData', reshape(xi_GaussZ(:,:,id),nlat*nlon,1));
        set(q, 'XData', xi_GaussX(:,:,id), ...
               'YData', xi_GaussY(:,:,id), ...
               'ZData', xi_GaussZ(:,:,id), ...
               'UData', u_GaussX(:,:,id), ...
               'VData', u_GaussY(:,:,id), ...
               'WData', u_GaussZ(:,:,id), 'Color', ColorInd);
        set(h, 'XData', CellX(:,:,id), ...
               'YData', CellY(:,:,id), ...
               'ZData', CellZ(:,:,id), ...
               'FaceAlpha', TransparencyInd);%, 'linestyle', 'none')
        set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
        axis(MinMaxBounds)

        %% Get frame as an image
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    close(v);
end