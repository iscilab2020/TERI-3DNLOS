function [rho_est, radiosity_tp, rho_hist, rad_hist, cost_fun] = EstRangeGivenAngleDistAccelerated(px,py,h,...
    Ntp,Ntp_sub,superFacets_ang,meas_ave,step_size,BGmodel,rho_0,rad_0, iter_max)
%FORWARDMODELDOORWAYCAMRANGE computes the contribution of light sources in the visible
% region (or visible side) to the visible (or relay) surface. This can be
% used a an approximate model for cancelling out background contributions.
%
%   Usage:
%       [h, Dh, rho_vals] = EstRangeGivenAngleDist(px,py,h,Ntp,Ntp_sub,superFacets_ang,rho_min_max)
%
%   Input:
%       * px (Mx-vector):	The measurement plane x-positions of the pixel locations.
%       * py (My-vector):	The measurement plane y-positions of the pixel locations.
%       * h  (scalar):  	The measurement plane x-positions of the pixel locations.
%       * Ntp (Mx-vector):	The measurement plane x-positions of the pixel locations.
%       * Ntp_sub (My-vector):	The measurement plane y-positions of the pixel locations.
%       * rho_hid  (scalar):  	The measurement plane x-positions of the pixel locations.
%
%   Output:
%       * Amat (MxMy-by-ns):        Matrix whose columns are contribution from the
%                                   n-th source.
%       * angles_tp (MxMy-by-ns):   Matrix whose columns are contribution from the
%                                   n-th source.
%       * angles_tp_sub (MxMy-by-ns):   Matrix whose columns are contribution from the
%                                       n-th source.
%
% Last Modified by $John Murray-Bruce$ at the University of South Florida
% v1.0 26-Feb-2023 15:17:32

max_iter = 1000; %500

if nargin<=9
    rho_0 = 0.5; 1.25;0.75; 0.5; 0.75;1.2;
    rho_0 = 0.75;
end
if nargin<=10
    rad_0 = 1;
end

% Set to true to see plot and values of range estimates with each iteration
diagnosticsConvergence = true;

superFacets_ang_bin = (superFacets_ang>1e-10);
[~, num_SF] = size(superFacets_ang_bin);

if nargin==7
    step_size = 5e-5;
elseif nargin == 8
    implicitBGsub = false;
elseif nargin == 9
    implicitBGsub = true;
    step_size_bg = 1e-8;
    bg_intensity = 0.01*[1 1]';
    BGmodel = BGmodel/max(BGmodel(:));
elseif nargin == 10
    implicitBGsub = true;
    step_size_bg = 1e-8;
    bg_intensity = 0.01*[1 1]';
    if isempty(BGmodel)
        implicitBGsub = false;
    end
elseif nargin==11
    implicitBGsub = true;
    step_size_bg = 1e-8;
    bg_intensity = 0.01*[1 1]';
    if isempty(BGmodel)
        implicitBGsub = false;
    end
elseif nargin==12
    implicitBGsub = true;
    step_size_bg = 1e-8;
    bg_intensity = 0.01*[1 1]';
    if isempty(BGmodel)
        implicitBGsub = false;
    end
    max_iter = iter_max;
end

if numel(step_size)==1
    step_size = [step_size step_size];
end

if isempty(rho_0)
    rho_0 = 0.75;
end

% Computational FOV
camFOV_edge_x = min(px);
min_angle_psi = atan(h/camFOV_edge_x);

Nt = Ntp(1);
Np = Ntp(2);
% Hidden volume discretization
psi_hid = linspace(min_angle_psi, pi/2-0.001, Np*Ntp_sub(2));
theta_hid = linspace(0.0,pi/2-0.001,Nt*Ntp_sub(1));

[TH, PH] = meshgrid(theta_hid,psi_hid);
TH = TH(:); PH = PH(:);
angles_tp_sub = [TH, PH];
superFacets_ang_mask = kron(reshape(sum(superFacets_ang_bin,2),Ntp),ones(Ntp_sub));
superFacets_ang_mask = superFacets_ang_mask(:);
lenSFmask = length(superFacets_ang_mask);
indexSF = (1:lenSFmask)';
indexSF = indexSF(superFacets_ang_mask>0);

if numel(rho_0)==1
    rho_t = rho_0*ones(1,num_SF); % 0.5 for others (1.5 only for Rf-Bc; Testing 0.1)
else
    rho_t = rho_0';
end

if numel(rad_0)==1
    rad_t = rad_0*ones(1,num_SF);
else
    rad_t = rad_0';
end
temp_vec = 1:prod(Ntp);
subFacet_indices = temp_vec(sum(superFacets_ang_bin,2)>=1);

for ii=1:num_SF
    rho_mask_temp = kron(reshape(superFacets_ang_bin(:,ii),Ntp),ones(Ntp_sub));
    rho_mask_tp(:,ii) = rho_mask_temp(:);
end
rhoVals_tp_mask = rho_mask_tp(indexSF,:);

% radiosity_tp = ones(1,num_SF);
radiosity_tp = rad_t';

radiosity_tp_prev = radiosity_tp;
rho_est_prev = rho_t';


rho_hist = zeros(num_SF,max_iter);
rad_hist = rho_hist;
% cost = zeros(1,max_iter);
y_rho = rho_est_prev;
u_rad = radiosity_tp_prev;
tp = 1;
superFacets_ang = superFacets_ang(subFacet_indices,:);
for ii = 1:max_iter
    rhoVals_tp = sum(rhoVals_tp_mask.*y_rho',2);
    [Am_rho_j, DelAm_rho_j] = fmDWRangeGivenAngleDist(px,py,h,subFacet_indices,angles_tp_sub,indexSF,Ntp,Ntp_sub,rhoVals_tp');

    temp_g_rad = (Am_rho_j*superFacets_ang);
    grho = sum(temp_g_rad .* u_rad',2);
    Dgrho = (DelAm_rho_j*superFacets_ang.* u_rad');

    if implicitBGsub
        grho = grho + BGmodel*bg_intensity;
        bg_intensity_k = bg_intensity - 2*step_size_bg*BGmodel'*(grho-meas_ave);
        bg_intensity = bg_intensity_k;
    end

    rho_est = max(y_rho - 2*step_size(1)*Dgrho'*(grho - meas_ave),0);
    rad_est_ktp = max(u_rad - 2*step_size(2)*temp_g_rad'*(grho - meas_ave),0);
%     rho_est = abs(y_rho - 2*step_size(1)*Dgrho'*(grho - meas_ave));
%     rad_est_ktp = abs(u_rad - 2*step_size(2)*temp_g_rad'*(grho - meas_ave));

    t = (sqrt(1+4*tp^2)+1)/2;
    %     alpha_k = (ii-1)./(ii+2);
    y_rho = rho_est + ((tp-1)/t).*(rho_est - rho_est_prev);
    u_rad = rad_est_ktp + ((tp-1)/t).*(rad_est_ktp - radiosity_tp_prev);
    rho_est_prev = rho_est;
    radiosity_tp_prev = rad_est_ktp;

    rho_hist(:,ii) = rho_est;
    rad_hist(:,ii) = rad_est_ktp;
    cost_fun(ii) = norm(grho-meas_ave)^2;

    rho_t = rho_est';
    radiosity_tp = rad_est_ktp';
    
    if (diagnosticsConvergence)&&(mod(ii,20)==0)
        ii
        figure(100); clf
        subplot(311)
        plot([1:ii],rho_hist(:,[1:ii])'); title('$\rho$ estimates')
        ylim([0 max([0.01;rho_hist(:)])]); xlim([0 max_iter])
        xlabel('iter')
        subplot(312)
        plot([1:ii],rad_hist(:,[1:ii])'); title('intensity estimates')
        xlabel('iter')
        ylim([0 max([0.01;rad_hist(:)])]); xlim([0 max_iter])
        subplot(313)
        plot([1:ii],cost_fun(1:ii));title('data fidelity')
        xlabel('iter')
        ylim([0 max([0.01;cost_fun(:)])]); xlim([0 max_iter])
        drawnow;
        rho_est

    end
end

end


function [Amat, DelAmat] = fmDWRangeGivenAngleDist(px,py,h,subFacet_indeces,angles_tp_sub,indexSF,Ntp,Ntp_sub,rhoVals_tp)
%FMDWRANGEGIVENANGLEDIST computes the measurements for given superfacets
%given its angular extent, assuming some fixed ranges rhoVals_tp dependent
%on the angles (theta and psi).

Mx = length(px);
My = length(py);
pz = h;
[Px, Py] = meshgrid(px,py);
meas_p = [Px(:), Py(:), pz*ones(Mx*My,1)];

TH = angles_tp_sub(indexSF,1)'; PH = angles_tp_sub(indexSF,2)';
denom_t = sqrt(1 + tan(TH).^2 + tan(PH).^2);
jacobian_factor = (sec(PH).^2.*sec(TH).^2)./(denom_t.^3);

[visMat] = dwvisibility_tp([TH' PH'], meas_p);
temp_VisJacobian = jacobian_factor.*visMat;

radial_falloff = ((rhoVals_tp./denom_t - meas_p(:,1)).^2 + ...
    (rhoVals_tp./denom_t.*tan(TH) - meas_p(:,2)).^2 + ...
    (rhoVals_tp.*tan(PH)./denom_t - meas_p(:,3)).^2 );
temp_rad = ((rhoVals_tp./denom_t - meas_p(:,1))./denom_t + ...
    (rhoVals_tp./denom_t.*tan(TH) - meas_p(:,2)).*(tan(TH)./denom_t ) + ...
    (rhoVals_tp.*tan(PH)./denom_t - meas_p(:,3)).*(tan(PH)./denom_t ) );

Del_radialfalloff = ((2*rhoVals_tp).*radial_falloff - (2*rhoVals_tp.^2).*(temp_rad))./(radial_falloff.^2);

Amat = ((rhoVals_tp.^2)./radial_falloff).*temp_VisJacobian;

DelAmat = Del_radialfalloff.*temp_VisJacobian;


% Nt = Ntp(1);
% Np = Ntp(2);
if prod(Ntp_sub)>1
    reduceMat_temp = zeros(prod(Ntp.*Ntp_sub),length(subFacet_indeces));
    for ii=1:length(subFacet_indeces)
        [subFacet_r, subFacet_c] = ind2sub(Ntp, subFacet_indeces(ii));    % subfacet 2D location (row, col)
        subFacet_rows_cols = ([subFacet_r, subFacet_c]-1).*Ntp_sub + (1:Ntp_sub)';
        [xa, xb] = meshgrid(subFacet_rows_cols(:,2), subFacet_rows_cols(:,1));
        reduceMat_temp(sub2ind(Ntp.*Ntp_sub, xb(:), xa(:)),ii) = 1;
    end
    reduceMat = reduceMat_temp(sum(reduceMat_temp,2)>=1,:);
    %
    %     reduceMat = kron(eye(Nt), kron(ones(Ntp_sub(1),1), kron(eye(Np), ones(Ntp_sub(2),1))));
    %     reduceMat = reduceMat(indexSF,:);
    Amat = Amat*reduceMat./prod(Ntp_sub);
    DelAmat = DelAmat*reduceMat./prod(Ntp_sub);
end

end


function [visMat] = dwvisibility_tp(s_tp,measurementplane)
%DWVISIBILITY computes and returns the visibility matrix of the point
% defined by s_rtp (in XZ-projected spherical coordinates) for all points
% in the measurement/observation plane.
%
%   Usage:
%       [visMat] = dwvisibility(s_rtp,measurementplane)
%
%   Input:
%       * s_rtp (TP-by-2 matrix):	    The spatial location of the hidden scene point
%                                       in XZ-projected spherical coordinates.
%       * measPlane (MxMy-by-3 matrix):	The xyz-locations on the measurement plane.
%
%   Output:
%       * visMat (MxMy-by-TP matrix):   Binary matrix indicating whether a hidden scene
%                                       point s_rtp is unoccluded (1) or occluded by the
%                                       doorway whose head is at z = 0, and coincides
%                                       with the y-axis.

% Last Modified by $John Murray-Bruce$ at the University of South Florida
% v1.0 26-Feb-2023 18:03:32

psi_s = s_tp(:,2)';
theta_s = s_tp(:,1)';
visibility_v = (( atan(measurementplane(:,2)./measurementplane(:,1)) - (theta_s)<0));
visibility_h = (( atan(measurementplane(:,3)./measurementplane(:,1)) - ( psi_s )<0));
visMat = visibility_h.*visibility_v;
end