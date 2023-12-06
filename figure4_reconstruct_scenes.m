%% Two-Edge-Resolved 3D NLOS Imaging (Figure 4 Reconstructions)
% Manuscript:
%   Czajkowski, R. and Murray-Bruce, J., 'Two-edge-resolved three-dimensional non-line-of-sight
% imaging with an ordinary camera', Nat. Commun., 2023.
% Notes:
% 1) Each newly computed forward model is saved, if it doesn't already
% exist. If you want to turn this off, set cacheForwardModels to false.

% Last Modified by John Murray-Bruce at the University of South Florida
% 03-Dec-2023 (Code clean-up and commented for dissemination)


%% Set configuration parameters and simulate forward model

clear; clc;
close all;

addpath('Functions')

cacheForwardModels = true;
% Function to convert from PESC to Cartesian coordinates for visualization
projsph2cartesian = @(rhoval, angles_tp) ([rhoval./sqrt(1 + tan(angles_tp(:,1)).^2 + tan(angles_tp(:,2)).^2), ...
    tan(angles_tp(:,1)).* rhoval./sqrt(1 + tan(angles_tp(:,1)).^2 + tan(angles_tp(:,2)).^2), ...
    tan(angles_tp(:,2)).* rhoval./sqrt(1 + tan(angles_tp(:,1)).^2 + tan(angles_tp(:,2)).^2)]);

% Height of doorway head
h = -(0.46*0.61);

% Number of measurement pixels (along each axis)
Mx = 125;
My = 125;

% Camera Field of View (small shift to avoid division by zero).
px = linspace(-0.46,0,Mx)-0.0000001;
py = linspace(-0.46,0,My)-0.0000001;

% Number of sub-discretizations used to simulate a each surface element.
Ntp_sub = [3 3];

N_theta = 120;      % Number of angular bins along azimuth (N_theta).
N_psi = 120;        % Number of angular bins along azimuth (N_theta).

rho_0 = 0.5;        % Initial range used for estimating shape.


fprintf('Simulating forward model...\n')
fwdmodelName = ['M',num2str(Mx),'x', num2str(My),'_N', num2str(N_theta), 'x',num2str(N_psi)];
fwdmodelpathName = ['ForwardModels', filesep, fwdmodelName,filesep, fwdmodelName,'rho_',num2str(rho_0),'.mat'];
% Simulate the forward model (or load it if it already exists).
tic
if isfile(fwdmodelpathName)
    load(fwdmodelpathName);
else
    [Amat,angles_tp] = forwardmodelDoorwayCam1(px,py,h,[N_theta N_psi], Ntp_sub,rho_0);
    if cacheForwardModels==true
        if ~isfolder(['ForwardModels', filesep, fwdmodelName,filesep])
            mkdir(['ForwardModels', filesep, fwdmodelName,filesep]);
        end
        save(fwdmodelpathName);
    end
end
toc

bgsource_loc = [-0.4 -0.4 2];
bg_term = simbackground(px,py,h,bgsource_loc);
bg_cols = [mean(Amat(:))*ones(My*Mx, 1) mean(Amat(:))*bg_term];
Bmat = [bg_cols Amat];

% Point to dataset path
addpath(genpath('Dataset'))
scene_names = {'','','','twoman0bg', 'twoman1bg', 'twoman2bg', 'twoman3bg'};

useTVrefinedRadiosity = true;


%% Manuscript Figure 4(a): Two Mannequins scene without visible side ambient light

scene_id = 4;           % Select 1 for play scene, 2 for work scene, and 3 for USF scene
fig_start = 0;

scriptRunReconstructionMethod

Ranges_2man0bg_scene = rho_est_t;

%% Manuscript Figure 4(b): Two Mannequins scene without visible side ambient light

scene_id = 5;           % Select 1 for play scene, 2 for work scene, and 3 for USF scene
fig_start = 4;

scriptRunReconstructionMethod

Ranges_2man1bg_scene = rho_est_t;

%% Manuscript Figure 4(c): Two Mannequins scene without visible side ambient light

scene_id = 6;           % Select 1 for play scene, 2 for work scene, and 3 for USF scene
fig_start = 8;

scriptRunReconstructionMethod

Ranges_2man2bg_scene = rho_est_t;

%% Display the range estimates
disp('======== RANGE ESTIMATES - TWO MANNEQUINS LOW BG  ========')
disp(['White marching mannequin: ', num2str(Ranges_2man0bg_scene(1)), 'm'])
disp(['Red-Green T-pose mannequin: ' num2str(Ranges_2man0bg_scene(2)),' m'])
fprintf('\n')

disp('======== RANGE ESTIMATES - TWO MANNEQUINS MODERATE BG  ========')
disp(['White marching mannequin: ', num2str(Ranges_2man1bg_scene(1)), 'm'])
disp(['Red-Green T-pose mannequin: ' num2str(Ranges_2man1bg_scene(2)),' m'])
fprintf('\n')

disp('======== RANGE ESTIMATES - TWO MANNEQUINS HIGH BG  ========')
disp(['White marching mannequin: ', num2str(Ranges_2man2bg_scene(3)), 'm'])
disp(['Red-Green T-pose mannequin: ' num2str(Ranges_2man2bg_scene(3)),' m'])
fprintf('\n')
