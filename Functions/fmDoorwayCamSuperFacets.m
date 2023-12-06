function [Amat, angles_tp, angles_tp_sub] = fmDoorwayCamSuperFacets(px,py,h,Ntp,Ntp_sub,rho_map)
%FMDOORWAYCAMSUPERFACETS computes the contribution of identified superfacets with 
% each superfacet identified with its estimated range from the corner.
%
%   Usage:
%       [Amat, angles_tp, angles_tp_sub] = forwardmodelDoorwayCam1(px,py,h,Ntp,Ntp_sub,rho_hid)
%       [Amat, angles_tp, angles_tp_sub] = forwardmodelDoorwayCam1(px,py,h,Ntp,Ntp_sub)
%       [Amat, angles_tp, angles_tp_sub] = forwardmodelDoorwayCam1(px,py,h,Ntp)
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



Nt = Ntp(1);
Np = Ntp(2);
if nargin==4
    Ntp_sub = [1 1];
end

Mx = length(px);
My = length(py);


rho_map_sub = kron(rho_map,ones(Ntp_sub));
Amat = zeros(Mx*My, Nt*Np*prod(Ntp_sub));
pz = h;
[Px,Py] = meshgrid(px,py);
meas_p = [Px(:),Py(:), pz*ones(Mx*My,1)];


% Computational FOV
camFOV_edge_x = min(px);
min_angle_psi = atan(h/camFOV_edge_x);

% Hidden volume discretization
% rho_hid = 1;
psi_hid = linspace(min_angle_psi, pi/2-0.001, Np*Ntp_sub(2));
theta_hid = linspace(0.0,pi/2-0.001,Nt*Ntp_sub(1));

[TH, PH] = meshgrid(theta_hid,psi_hid);
TH = TH(:); PH = PH(:);
angles_tp_sub = [TH, PH];

denom_t = sqrt(1 + tan(TH).^2 + tan(PH).^2);
jacobian_factor = (sec(PH).^2.*sec(TH).^2)./(denom_t.^3);

for ii = 1:length(denom_t)
    rho_hid = rho_map_sub(ii);
    radial_falloff = (rho_hid.^2)./( (rho_hid./denom_t(ii) - meas_p(:,1)).^2 + ...
        (rho_hid./denom_t(ii).*tan(TH(ii)) - meas_p(:,2)).^2 + ...
        (rho_hid./denom_t(ii).*tan(PH(ii)) - meas_p(:,3)).^2 );
    [visMat] = dwvisibility([1 angles_tp_sub(ii,:)],meas_p);
    Amat(:,ii) = radial_falloff.*jacobian_factor(ii).*visMat;

end

if prod(Ntp_sub)>1
    reduceMat = kron(eye(Nt), kron(ones(Ntp_sub(1),1), kron(eye(Np), ones(Ntp_sub(2),1))));
    Amat = Amat*reduceMat./prod(Ntp_sub);
end

[TH, PH] = meshgrid(linspace(0.000,pi/2-0.001,Nt),...
    linspace(min_angle_psi, pi/2-0.001,Np));
angles_tp = [TH(:),PH(:)];


end



function [visMat] = dwvisibility(s_rtp,measurementplane)
%DWVISIBILITY computes and returns the visibility matrix of the point
% defined by s_rtp (in XZ-projected spherical coordinates) for all points
% in the measurement/observation plane.
%
%   Usage:
%       [visMat] = dwvisibility(s_rtp,measurementplane)
%
%   Input:
%       * s_rtp (3-vector):	            The spatial location of the hidden scene point 
%                                       in XY-projected spherical coordinates. 
%       * measPlane (MxMy-by-3 matrix):	The xyz-locations on the measurement plane.
%
%   Output:
%       * visMat (MxMy-vector):         Binary vector indicating whether a hidden scene
%                                       point s_rtp is unoccluded (1) or occluded by the 
%                                       doorway whose head is at z = 0, and coincides
%                                       with the y-axis.

% Last Modified by $John Murray-Bruce$ at the University of South Florida
% v1.0 26-Feb-2023 18:03:32

psi_s = s_rtp(3);
theta_s = s_rtp(2);
rho_s = s_rtp(1);

x = rho_s./sqrt(1 + tan(theta_s).^2 + tan(psi_s).^2);
y = x * tan(theta_s);
z = x*tan(psi_s);

s = [x y z];


visibility_v = (( atan(measurementplane(:,2)./measurementplane(:,1)) - (theta_s)<0));
visibility_h = (( atan(measurementplane(:,3)./measurementplane(:,1)) - ( psi_s )<0));

visMat = visibility_h.*visibility_v;

end