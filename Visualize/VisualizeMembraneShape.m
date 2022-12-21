%% Visualize membrane shape while flowing in an unbounded domain
%% Create grid for representation: [N+1 x 2*N+1]
nlat = N+1;
nlon = 2*N+1;

%% Visualization settings 
VisualizeSettings

%% Patch faces (connectivity matrix)
faces = PatchFaces(nlat, nlon);

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
Vertices = zeros((nlat+2)*nlon,3,Numframe);

nframe = 0;
T_step = zeros(Numframe,1);
for nstep = 0:NSTEPS-1
    %% Read from file
    cxi = fread(fidCoord,3*(N+1)^2,'double');
    axi = zeros(size(mask_a));  bxi = zeros(size(mask_b));
    axi(mask_a) = cxi(1:3*(N+1)*(N+2)/2);
    bxi(mask_b) = cxi(3*(N+1)*(N+2)/2+1:3*(N+1)^2);

    if (nstep~=0 && nstep~=NSTEPS-1 && mod(nstep,timeStepIncrement)~=0)
        continue
    end
    nframe = nframe + 1;
    T_step(nframe) = Time(nstep+1);
    
    %% Coordinates of vertices
    xi_equal = shsecm(axi,bxi);
    xi_Gauss = shsgcm(axi,bxi);
    xiplot = [xi_equal(1,:,:); xi_Gauss; xi_equal(end,:,:)];
    Vertices(:,:,nframe) = reshape(xiplot,(nlat+2)*nlon,3);
end
%% Dimensionalization 
T_step = T_step/RefShearRate; % in seconds
Vertices = Vertices*RefLength*10^(6); % \mum
%%
MinMaxBounds = [min(min(Vertices(:,1,:))) max(max(Vertices(:,1,:))) ...
                min(min(Vertices(:,2,:))) max(max(Vertices(:,2,:))) ...
                min(min(Vertices(:,3,:))) max(max(Vertices(:,3,:)))];
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

hRBC = patch('Faces',faces','Vertices',Vertices(:,:,1),...
             'FaceColor','r','FaceAlpha',0.6,...
             'EdgeColor',ColorInd, 'EdgeAlpha',1,'FaceLighting','gouraud',...
             'SpecularColorReflectance',0.1);
set(gca,'FontName','cambria math','FontSize',12)
camlight
material Dull
set(gca,'DataAspectRatio',[1 1 1])
ht = title({sprintf('Time = %4.4f sec',T_step(1))},'Color',ColorInd);
set(ht,'FontName','cambria math','FontSize',12)
view(viewInd)
axis(MinMaxBounds)
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
        set(hRBC, 'Vertices', Vertices(:,:,id));
        set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
        axis(MinMaxBounds)
        
        %% Get frame as an image
        frame = getframe(gcf);

        %% Create a colormap for the first frame. For the rest of the frames,
        %% use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(frame.cdata, 1024, 'dither');
        else
            mov(:,:,1,id) = rgb2ind(frame.cdata, map, 'dither');
        end
    end
    %% Create animated GIF
    imwrite(mov, map, 'MembraneShape.gif', 'DelayTime', 0, 'LoopCount', inf)
else % Write to MP4 file 
    v = VideoWriter('MembraneShape.mp4','MPEG-4');
    v.FrameRate = 10;
    open(v);

    %% Loop through by changing XData and YData
    for id = 1:length(T_step)
        %% Update graphics data. This is more efficient than recreating plots.
        set(hRBC, 'Vertices', Vertices(:,:,id));
        set(ht, 'String', {sprintf('Time = %4.4f sec',T_step(id))})
        axis(MinMaxBounds)

        %% Get frame as an image
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    close(v);
end