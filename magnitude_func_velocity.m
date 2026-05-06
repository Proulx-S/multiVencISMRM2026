function m = magnitude_func_velocity(v, B, C1, C2, C3, C4)
% Polynomial magnitude as a function of velocity: m(v) = B + C1*v + C2*v^2 + ...
% C3 and C4 are optional (default 0).
if nargin < 5; C3 = 0; end
if nargin < 6; C4 = 0; end
m = max(0,   B + C1.*v + C2.*v.^2 + C3.*v.^3 + C4.*v.^4   );
