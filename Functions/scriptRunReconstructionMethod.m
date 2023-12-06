%% Two-edge-resolved NLOS imaging
% ------------------------------Pseudo-code------------------------------
% STEP 1. Shape-only reconstruction with L1/sparsity prior.
%           (a) A linear inverse problem (LIP) is solved, per colour channel,
%                to recover the azimuth and projected-elevation representation 
%                (i.e., the shapes) of hidden scene objects, with hidden scene
%                assumed to be confined to a single fixed range.
%           (b) The reconstruction is analysed for connected surface elements
%                that most likely belong to the same cluster. 
%
% STEP 2. Range estimation and scene radiosity refinement with smoothness prior.
%           (c) A non-linear inverse problem (NLIP) is solved to estimate the
%                ranges and global radiosities of each identified cluster
%                of connected elements.
%           (d) A final 3D colour reconstruction is produced by incorporating
%                the estimated ranges and solving a resulting total-variation
%                -regularized LIP.
% -----------------------------------------------------------------------

% Manuscript:
%   Czajkowski, R. and Murray-Bruce, J., 'Two-edge-resolved three-dimensional non-line-of-sight
% imaging with an ordinary camera', Nat. Commun., 2023.

% Last Modified by John Murray-Bruce at the University of South Florida
% 03-Dec-2023 (Code clean-up and commented for dissemination)

% To do list:
% 1) use 'filesep' to make code compatible with Windows and Linux
% 2) Add option to save or not save forward model


% Set parameters for each scene configuration!
if scene_id == 1
    % PLAY scene

    % Set l1-reg parameter
    lambda_val_l1 = 25/100; % For Ntp = [120 120]; --- Play Scene


    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1500;

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010;

    % Set gamma
    thresh_val_1 = 0.003;

    % Set J_max
    maxNumObjs = 4;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = .945*[6 10.5]*(1e-5);
    rho_start = rho_0+0.5;
    rad_start = 1; % rad_est'/(max(rad_est));
    num_iterations = 2000;

    rho_map = 5*ones(N_theta,N_psi);

    try
        if singleSnapshotData==true
            step_size = .5*[11 8]*(1e-5);
        end
    catch
    end
elseif scene_id == 2
    % WORK scene

    % Set l1-reg parameter
    lambda_val_l1 = 20/100; % WORK --- For Ntp = [120 120]

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1500;

    lambda_val_tv = .00000000000010/5;

    % Set gamma
    thresh_val_1 = 0.000;

    % Set J_max
    maxNumObjs = 4;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = 9.3*[.09*(1e-5) .035*(1e-5)]; % WORK - Good
    rho_start = rho_0+1;
    rad_start = 1;
    num_iterations = 500;

    rho_map = 5*ones(N_theta,N_psi);
elseif scene_id == 3
    % USF+Man scene

    % Set l1-reg parameter
    lambda_val_l1 = 46/100; % USF --- For Ntp = [120 120]

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1500;

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010;

    % Set gamma
    thresh_val_1 = 0.002;

    % Set J_max
    maxNumObjs = 5;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = 1*[.09*(1e-3)/100 .035*(1e-3)/100]; % USF
    rho_start = rho_0+0.75;
    rad_start = 1;
    num_iterations = 500;

    rho_map = 5*ones(N_theta,N_psi);
elseif scene_id == 4
    % Two Man scene without visible side ambient illumination

    % Set l1-reg parameter
    lambda_val_l1 = 1/15*[4 3 4];

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1000;

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010/5;

    % Set gamma
    thresh_val_1 = 0.00;

    % Set J_max
    maxNumObjs = 3;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = 3*[.9*(1e-5) .35*(1e-5)]; % 2Man0BG
    rho_start = rho_0+0.5;
    rad_start = 1;
    num_iterations = 500;

    rho_map = 10*ones(N_theta,N_psi);
elseif scene_id == 5
    % Two Man scene with moderate amount of visible side ambient illumination

    % Set l1-reg parameter
    lambda_val_l1 = 1/30*[4 3 4];

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1000;  % Two Man 1BG

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010/5;

    % Set gamma
    thresh_val_1 = 0.00;

    % Set J_max
    maxNumObjs = 3;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = 1*[.9*(1e-5) .35*(1e-5)]; % 2Man1BG
    rho_start = rho_0+0.5;
    rad_start = 1;
    num_iterations = 500;

    scale_color = 50;
    rho_map = 10*ones(N_theta,N_psi);
elseif scene_id == 6
    % Two Man scene with high amount of visible side ambient illumination

    % Set l1-reg parameter
    lambda_val_l1 = 1*(1/2000)*[4 3 3]; % Two Man BG2

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1000;

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010/5;

    % Set gamma
    thresh_val_1 = 0.000;

    % Set J_max
    maxNumObjs = 3;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = .5*[.9*(1e-5) 1*(1e-5)];
    rho_start = rho_0+0.5;
    rad_start = 1;
    num_iterations = 500;
    rho_map = 10*ones(N_theta,N_psi);
elseif scene_id == 7
    % Two Man scene with high amount of visible side ambient illumination

    % Set l1-reg parameter
    lambda_val_l1 = 0.5*(1/2000)*[4 3 3];

    % Set maximum number of iterations L1-reg
    MAX_it_l1 = 1000;  % Two Man BG2

    % Set tv-reg parameter
    lambda_val_tv = .00000000000010/5;

    % Set gamma
    thresh_val_1 = 0.00;

    % Set J_max
    maxNumObjs = 3;

    % Set gradient algorithm parameters for range estimation via NLIP
    % Set step_size
    step_size = .5*[.9*(1e-5) .2*(1e-5)]; % 2Man2BG
    rho_start = rho_0+0.5;
    rad_start = 1;
    num_iterations = 500;
    rho_map = 10*ones(N_theta,N_psi);
else
    fprintf('Invalid scene ID\n')
end


%% Load measured photograph and perform shape-only reconstruction (STEP 1)

disp('Step 1: (Fig. 2a) Loading penumbra photograph')
% Load red channel
im_r = load([scene_names{scene_id}, '_red.mat']);
im_r = flipud(flipud(im_r.im)');
real_meas_r = im_r(:);

% Load green channel
im_g = load([scene_names{scene_id},'_green.mat']);
im_g = flipud(flipud(im_g.im)');
real_meas_g = im_g(:);

% Load blue channel
im_b = load([scene_names{scene_id},'_blue.mat']);
im_b = flipud(flipud(im_b.im)');
real_meas_b = im_b(:);
figure(fig_start+1);
image(cat(3,rot90(flipud(im_r)',2),rot90(flipud(im_g)',2),rot90(flipud(im_b)',2)));
title('Penumbra photograph','FontSize',20)
disp('Step 1: (Fig. 2b) Commencing shape reconstruction per color channel')
tic
% Perform shape-only reconstruction with range initialized to 0.5.
auglag_par = 0.99;
overrelax_par = 1.05;
[z, ~] = lassoADMM_RGB(Bmat, real_meas_r, real_meas_g, real_meas_b, lambda_val_l1, ...
    auglag_par, overrelax_par, MAX_it_l1);
scene_est_l1_bg1_r = z.r;
scene_est_l1_bg1_g = z.g;
scene_est_l1_bg1_b = z.b;

disp(['      *Shape reconstruction time taken [s]: ',num2str(toc)])

scene_est_l1_bg_r = scene_est_l1_bg1_r(3:end);
scene_est_l1_bg_g = scene_est_l1_bg1_g(3:end);
scene_est_l1_bg_b = scene_est_l1_bg1_b(3:end);

bgContribution = bg_cols;

%% Set gamma parameter, i.e. threshold value.

% Plot angle-estimation solution
% Plotting parameters
L_max = 0.65;               % Max limit of coordinate axes.
font_size = 20;             % Font size of plot texts.

% Form binary mask that will be used for clustering objects
scene_est_l1_bg_binMask = zeros(N_theta*N_psi,1);
scene_est_l1_bg_binMask(scene_est_l1_bg_r>thresh_val_1,1) = 1;
scene_est_l1_bg_binMask(scene_est_l1_bg_g>thresh_val_1,1) = 1;
scene_est_l1_bg_binMask(scene_est_l1_bg_b>thresh_val_1,1) = 1;
est_ang_pix = (angles_tp(scene_est_l1_bg_binMask>thresh_val_1,:));
scene_est_xyz = projsph2cartesian(rho_0,est_ang_pix);
scene_est_temp = [scene_est_l1_bg_r scene_est_l1_bg_g scene_est_l1_bg_b];
scene_est_temp(scene_est_temp<thresh_val_1) = 0;


disp('Step 1: (Fig. 2c) Plotting the shape-only reconstruction in colour')
% Plot shape-only reconstruction in colour (scaling '5/50' used here to enhance
% visualization, not used in in reconstruction in any way).
sceneShape_col = scene_est_temp(scene_est_l1_bg_binMask>thresh_val_1,:)./max(scene_est_temp(:))*5;
if scene_id>=5
    sceneShape_col = scene_est_temp(scene_est_l1_bg_binMask>thresh_val_1,:)./max(scene_est_temp(:))*40;
end
figure(fig_start+2); clf;
scatter3(scene_est_xyz(:,1),scene_est_xyz(:,2),scene_est_xyz(:,3),10,sceneShape_col,'filled');
hold on;
axis equal
xlim([0 L_max])
ylim([0 L_max])
zlim([0 L_max])
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)
xlabel('$x\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
ylabel('$y\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
zlabel('$z\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
grid on;
grid minor
box on;


% This shows the different colour channels individually (as in manuscript Fig. 2b)
plotColChannels = false;
if plotColChannels == true
    script_PlotEstimates_perColourChannel;
end
disp('Step 1: (Fig. 2d) Identifying all clusters in the shape-only reconstruction')

% Perform clustering into clusters of surface elements
[superFacets,CC,superFacets_bin] = computeSuperFacetsOfSceneEstimate(...
    scene_est_l1_bg_binMask,[N_theta N_psi],0,maxNumObjs,true);

% Plot the indentified clusters.
n_subplots = ceil(sqrt(maxNumObjs));
rad_est = [];
str_leg = {};
figure(fig_start+3);clf;
for ii=1:size(superFacets,2)
    rad_est(ii) = sum(sum(scene_est_temp(superFacets(:,ii)>0,:),2)/3)./sum(superFacets(:,ii));
    figure(fig_start+3);
    subplot(n_subplots,ceil(maxNumObjs/n_subplots),ii)
    imagesc(reshape(superFacets(:,ii),[N_theta N_psi]))
    str_leg{ii} = num2str(ii);
    axis equal square tight
    axis off
    colormap hot
end

rad_est_pix = sum(scene_est_temp,2);
if scene_id>3
    rad_est_pix = rad_est;
end
drawnow;



%% STEP 2: Range recovery
% Recover range by solving non-linear inverse problem

disp('Step 2: (Fig. 2e) Commencing gradient algorithm for per cluster range estimation')

if scene_id == 1
    % Play scene uses this initial radiosity in gradient algorithm (no other uses this)
    rad_start = rad_est'/(max(rad_est));
end
tic
% Run gradient algorithm to estimate ranges for each cluster
real_meas_ave = (real_meas_r + real_meas_g + real_meas_b)/3;
[rho_est_t, rad_est_t, rho_hist, rad_hist,cost_fun] = EstRangeGivenAngleDistAccelerated(...
    px, py, h, [N_theta N_psi], [1 1], superFacets.*(rad_est_pix),...
    real_meas_ave, step_size,[],rho_start,rad_start,num_iterations);
disp(['      *Range estimation time taken [s]: ',num2str(toc)])

%% Plot range and shape solution (without colour)
L_max = 1.0;
II = min(CC.NumObjects,size(superFacets,2));
figure(fig_start+3); clf;
if scene_id == 6
    % For this high ambient light scene, the last cluster is actually
    %  the true scene and not spurious clutter.
    % For all other scenes the last cluster is unwanted clutter, and so we
    % opt to not display them. In this case, we display everything.
    II = II+1;
end
for ii=1:II-1
    est_pix_xyz = projsph2cartesian(rho_est_t(ii), angles_tp(superFacets(:,ii)>0,:));
    Col = superFacets(superFacets(:,ii)>0,ii)*ones(1,3)./max(superFacets(:,ii))...
        *mean(rad_est_t(ii))*3;
    scatter3(est_pix_xyz(:,1),est_pix_xyz(:,2),est_pix_xyz(:,3),20,Col,'filled');
    hold on;
    if ~isnan(rho_est_t(ii))
        rho_map(superFacets(:,ii)>0) = rho_est_t(ii);
    end
end
box on
axis equal
xlim([0 L_max])
ylim([0 L_max])
zlim([0 L_max])
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)
xlabel('$x\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
ylabel('$y\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
zlabel('$z\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
grid on;
grid minor
box on;

%% STEP 2: Refine radiosity estimate via TV regularization
% Enforces smoothness prior (to obtain better colour renditions)

if useTVrefinedRadiosity == true
    % Main manuscript algorithm (uses radiosity refinement via
    % TV/smoothness prior

    disp('Step 2: (Fig. 2f) Commencing radiosity refinement with a TV prior')


    disp('    *Computing new FWD Model at new ranges')
    tic
    % Compute forward model
    if scene_id<4
        % Can use range estimate or a constant range to form the forward
        % model (with negligible impact: an effect of info orthogonality)
        [Arho_mat] = fmDoorwayCamSuperFacets(px,py,h,[N_theta N_psi],Ntp_sub,ones(size(rho_map)));
    else
        % This uses the estimated ranges to form the forward model
        [Arho_mat] = fmDoorwayCamSuperFacets(px,py,h,[N_theta N_psi],Ntp_sub,rho_map);
    end
    disp(['    *FWD Model computation time taken [s]: ',num2str(toc)])


    tt = tic;
    disp('    *Refining radiosity estimates')
    % TV reconstruction
    [reconstruction] = fistatvrecon3DCC(Arho_mat, [real_meas_r, real_meas_g, real_meas_b], ...
        lambda_val_tv, [N_theta N_psi]);
    
    % Gather solutions of TV reconstruction
    scene_est_tv_second_r = reconstruction(:,:,1);
    scene_est_tv_second_g = reconstruction(:,:,2);
    scene_est_tv_second_b = reconstruction(:,:,3);
    scene_est_tv_second_r = scene_est_tv_second_r(:);
    scene_est_tv_second_g = scene_est_tv_second_g(:);
    scene_est_tv_second_b = scene_est_tv_second_b(:);

    disp(['    *FISTA_TV radiosity refinement time taken [s]: ',num2str(toc(tt))])
end
%% Plot scene reconstruction/estimate
% Display the computed solution (global scalings and thresholds chosen
% ad hoc are used to enhance display of the estimated scenes)
if scene_id==3
    L_max = 1.2;
end
figure(fig_start+4); clf;


if useTVrefinedRadiosity==true
    scene_est_tv_ref_r = scene_est_tv_second_r;
    scene_est_tv_ref_g = scene_est_tv_second_g;
    scene_est_tv_ref_b = scene_est_tv_second_b;
    thresh_val = 0.000;%0.0002
    if scene_id==4
        % Tweaks to optimize display of scene
        thresh_val = 0.000125;
    elseif scene_id == 6
        thresh_val = 0.0000;
    end
    % Only interested in plotting points with radiosity larger than some threshold
    scene_est_tv_ref = zeros(N_theta*N_psi,1);
    scene_est_tv_ref(scene_est_tv_ref_r>thresh_val,1) = 1;
    scene_est_tv_ref(scene_est_tv_ref_g>thresh_val,1) = 1;
    scene_est_tv_ref(scene_est_tv_ref_b>thresh_val,1) = 1;
    est_ang_pix = (angles_tp(scene_est_tv_ref>thresh_val,:));
    rho_map_flat = rho_map(:);
    rho_map_flat = rho_map_flat(scene_est_tv_ref>thresh_val);
    scene_est_xyz = projsph2cartesian(rho_map_flat,est_ang_pix);
    scene_est_temp = [scene_est_tv_ref_r scene_est_tv_ref_g scene_est_tv_ref_b];
    scene_est_temp(scene_est_temp<thresh_val) = 0;
    Col = scene_est_temp(scene_est_tv_ref>thresh_val,:)./max(scene_est_temp(:))*20;

    if scene_id == 3
        % Tweaks to optimize display of scene
        % Can get better visualization if each object is scalled separately
        scene_est_temp(scene_est_temp<thresh_val) = 0;
        rad_temp = rad_est_t;
        rad_temp(3) = rad_temp(3)*0.5;
        rad_temp(4:5) = rad_temp(4:5)*3;
        scene_est_temp = scene_est_temp.*sum(superFacets.*rad_temp,2);
        Col = scene_est_temp(scene_est_tv_ref>thresh_val,:)./max(scene_est_temp(:))*14;
    elseif scene_id >4
        % Scaling to optimize display of scene
        Col = scene_est_temp(scene_est_tv_ref>thresh_val,:)./max(scene_est_temp(:))*80;
    end
    if scene_id ==6
        % Scaling to optimize display of scene
        Col = scene_est_temp(scene_est_tv_ref>thresh_val,:)./max(scene_est_temp(:))*200;
    end
else
    thresh_val = 0.00;
    scene_est_l1_bg_ref(scene_est_l1_bg_r>thresh_val,1) = 1;
    scene_est_l1_bg_ref(scene_est_l1_bg_g>thresh_val,1) = 1;
    scene_est_l1_bg_ref(scene_est_l1_bg_b>thresh_val,1) = 1;
    est_ang_pix = (angles_tp(scene_est_l1_bg_ref>thresh_val,:));
    rho_map_flat = rho_map(:);
    rho_map_flat = rho_map_flat(scene_est_l1_bg_ref>thresh_val);
    scene_est_xyz = projsph2cartesian(rho_map_flat,est_ang_pix);
    scene_est_temp = [scene_est_l1_bg_r scene_est_l1_bg_g scene_est_l1_bg_b];
    scene_est_temp(scene_est_temp<thresh_val) = 0;
    Col = scene_est_temp(scene_est_l1_bg_ref>thresh_val,:)./max(scene_est_temp(:))*5.5;
end

% Scatter3 is used to plot point cloud (could plot surface elements instead).
scatter3(scene_est_xyz(:,1),scene_est_xyz(:,2),scene_est_xyz(:,3),12*(.125+rho_map_flat).^2,Col,'filled');
hold on;
box on
axis equal
xlim([0 L_max+0/4])
ylim([0 L_max+0/4])
zlim([0 L_max+0/4])

set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)

xlabel('$x\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
ylabel('$y\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
zlabel('$z\, [{\rm m}]$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
grid on;
grid minor
box on;

