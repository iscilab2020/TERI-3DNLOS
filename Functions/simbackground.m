function [Bmat] = simbackground(px,py,h,locs_s)
%SIMBACKGROUND computes the contribution of light sources in the visible
% region (or visible side) to the visible (or relay) surface. This can be
% used a an approsimate model for cancelling out background contributions.

%   Usage:
%       [bmat] = simbackground(px, py, h, bgsource)
%
%   Input:
%       * px (Mx-vector):	The measurement plane x-positions of the pixel locations.
%       * py (My-vector):	The measurement plane y-positions of the pixel locations.
%       * h  (scalar):  	The measurement plane x-positions of the pixel locations.
%       * locs_s (ns-by-3):	Source locations in Cartesian coordinates (ns sources).
%
%   Output:
%       * Bmat (MxMy-by-ns):Matrix whose columns are contribution from the
%                           n-th source.

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 26-Feb-2023 15:17:32


Mx = length(px);
My = length(py);

[ns, ~] = size(locs_s);

pz = h;
[Px,Py] = meshgrid(px,py);
meas_p = [Px(:),Py(:), pz*ones(Mx*My,1)];

for ii = 1:ns
    Bmat(:,ii) = (1./sum((meas_p - ones(My*Mx,1)*locs_s(ii,:)).^2,2));
end