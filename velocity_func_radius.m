function v = velocity_func_radius(r, Vmax, R)
% Parabolic velocity profile: v(r) = Vmax * (1 - (r/R)^2), clamped to 0 at the wall.
v = max(0,   Vmax .* (1 - (r ./ R).^2)   );
