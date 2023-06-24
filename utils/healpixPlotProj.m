function healpixPlotProj(map,projType,lookingOut)
% function healpixPlotProj(map,projType,lookingOut)
% map is a healpix pixel map to plot
% projType is the projection type to plot
%   Supports:
%    - TrueView
% lookingOut is a logical to specify if we are looking from the centre of
% the sphere out (true), or outside the sphere in (false). Typically global 
% sky maps are view by looking out of the sphere - we are in the centre of
% the universe... So this is the default
% useEqArea uses the MAAT equal area version of the projection which is different from the standard one (false default)

if nargin < 2 || isempty(projType), projType = 'TrueView'; end
if nargin < 3, lookingOut = true; end
signPhi = 1;
if lookingOut, signPhi = -1; end

sz = size(map);
nside = sqrt(max(sz)/12);

% Use the MEALpix package
tp = pix2ang(nside);
tp = [tp{:}];
th = tp(1,:);
ph = signPhi.*tp(2,:);

[u,v,w] = PhTh2DirCos(ph,th);

switch projType
    case 'TrueView'
        [x,y] = DirCos2TrueView(u,v,w);
    otherwise
        error(['Unknown projection type: ', projType])
end

gridDelaunay = delaunay(x,y);
% h = trisurf(gridDelaunay,x,y,map*0.0,map);
h = trisurf(gridDelaunay,x,y,map,map);

set(h, 'LineStyle', 'none')
axis equal
axis off
view([0,90])


