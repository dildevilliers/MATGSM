function healpixPlotMollweide(map,lookingOut,useEqArea)
% function healpixPlotMollweide(map,lookingOut,useEqArea)
% map is a healpix pixel map to plot
% lookingOut is a logical to specify if we are looking from the centre of
% the sphere out (true), or outside the sphere in (false). Typically global 
% sky maps are view by looking out of the sphere - we are in the centre of
% the universe... So this is the default
% useEqArea uses the MAAT equal area version of the projection which is different from the standard one (false default)

if nargin < 3 || isempty(useEqArea), useEqArea = false; end

if nargin < 2, lookingOut = true; end
signPhi = 1;
if lookingOut, signPhi = -1; end

sz = size(map);
nside = sqrt(max(sz)/12);

% Use the MEALpix package
tp = pix2ang(nside);
tp = [tp{:}];
th = tp(1,:);
ph = signPhi.*tp(2,:);

if useEqArea
    [x, y] = PhTh2MollweideEqArea(ph, th);
else
    [x, y] = PhTh2Mollweide(ph, th);
end

gridDelaunay = delaunay(x,y);
% h = trisurf(gridDelaunay,x,y,map*0.0,map);
h = trisurf(gridDelaunay,x,y,map,map);

set(h, 'LineStyle', 'none')
axis equal
axis off
view([0,90])


