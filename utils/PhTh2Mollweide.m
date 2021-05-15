function [x, y] = PhTh2Mollweide(ph, th)
% [x, y] = PhTh2Mollweide(ph, th)

MAX_ITER = 1e5;
TOL = 1e-10;

% Convert theta to longitude.
th = pi/2 - th(:);
% Wrap ph to +-pi
ph = wrap2pi(ph(:));

t = th;
for it = 1:MAX_ITER
   dt = (t + sin(t) - pi.*sin(th))./(1 + cos(t));
   t = t - dt;
   if(max(abs(dt)) < TOL)
      break;
   end
end
t = t/2;
x = 2.*sqrt(2)./pi.*ph.*cos(t);
y = sqrt(2).*sin(t);