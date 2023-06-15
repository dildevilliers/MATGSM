function [x, y] = PhTh2MollweideEqArea(ph, th)
% [x, y] = PhTh2MollweideEqArea(ph, th)
% This below is the MAAT equal area Mollweide version

ph = wrap2pi(ph);
th = pi/2 - th;
x = 2.*sqrt(2)./pi.*ph.*cos(th);
y = sqrt(2).*sin(th);

