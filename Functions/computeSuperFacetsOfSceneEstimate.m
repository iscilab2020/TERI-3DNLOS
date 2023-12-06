function [superFacets,CC,superFacets_hull] = computeSuperFacetsOfSceneEstimate(scene_est,Ntp,th_val,maxNumObj,downsampMask)
%COMPUTESUPERFACETSOFSCENEESTIMATE combines contiguous facet/surface elements of the
% estimated scene and treats them as a single hidden scene element (or target).
% By contiguous, we mean angularly neighboring elements as computed using
% connected component analysis (i.e., using MATLAB's built-in bwconncomp()
% function).
%
%   Usage:
%       [superFacets,CC] = computeSuperFacetsOfSceneEstimate(scene_est,Ntp,th_val)
%
%   Input:
%       * scene_est (N-vector):	An estimate of the hidden scene's radiosity wrt angles.
%       * Ntp (2-vector):	    Hidden scene discretization over angles.
%       * th_val  (scalar):  	Threshold level for removing spurious hidden scene estimates.
%
%   Output:
%       * superFacets (MxMy-by-ns):        Matrix whose n-th column models the contribution of the
%                                   n-th hidden scene facet/patch.
%       * CC (cell):    The angular spatial locations of the
%                                   surface elements: i.e., (theta, psi).
%       * angles_tp_sub (MxMy-by-ns):   Matrix whose columns are contribution from the
%                                       n-th source.
%
% Last Modified by $John Murray-Bruce$ at the University of South Florida
% v1.0 26-Feb-2023 15:17:32

if nargin==3
    maxNumObj = 5;
    downsampMask = false;
elseif nargin == 4
    downsampMask = false;
end


binMask = scene_est>th_val;
binMask = reshape(binMask,Ntp);
if downsampMask==true
    binMask = conv2(binMask,0.5*[ 1 1; 1 1],'same');
end
CC = bwconncomp(binMask,26);
% imagesc(binMask);

imageSegments = zeros(size(binMask));
lenImageSegs = zeros(CC.NumObjects,1);
for ii=1:CC.NumObjects
    lenImageSegs(ii) = length(CC.PixelIdxList{ii});
    imageSegments(CC.PixelIdxList{ii}) = lenImageSegs(ii);
end

[~, indx] = sort(lenImageSegs,'descend');

maxNumObj = min(maxNumObj,CC.NumObjects);
superFacets = zeros(length(scene_est),maxNumObj);
superFacets_hull = zeros(length(scene_est),maxNumObj);
for ii=1:maxNumObj-1
    SF_bin_temp = zeros(Ntp);
    superFacets(CC.PixelIdxList{indx(ii)},ii) = scene_est(CC.PixelIdxList{indx(ii)},1);
    SF_bin_temp(CC.PixelIdxList{indx(ii)}) = 1;
    SF_bin_temp = bwconvhull(SF_bin_temp,'union');
    superFacets_hull(:,ii) = SF_bin_temp(:);
end

superFacets(vertcat(CC.PixelIdxList{indx(maxNumObj:end)}),maxNumObj) = ...
                            scene_est(vertcat(CC.PixelIdxList{indx(maxNumObj:end)}),1);

superFacets_hull(vertcat(CC.PixelIdxList{indx(maxNumObj:end)}),maxNumObj) = 1;

end