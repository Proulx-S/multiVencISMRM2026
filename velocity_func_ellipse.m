function v = velocity_func_ellipse(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset)
% Parabolic velocity profile with elliptical wall centred at velocity peak.
% Peak = ellipse centre at (FEoffset, PEoffset).
% Semi-major axis R along angle alpha from PE axis; semi-minor axis R/AR (AR>=1).
dPE = r.*cos(p)  - PEoffset;
dFE = -r.*sin(p) - FEoffset;
r_v = sqrt(dPE.^2 + dFE.^2);
uPE = dPE ./ max(r_v, eps);
uFE = dFE ./ max(r_v, eps);
A     = uPE.*cos(alpha) + uFE.*sin(alpha);
B     = -uPE.*sin(alpha) + uFE.*cos(alpha);
R_eff = R ./ sqrt(A.^2 + AR.^2.*B.^2);
v = Vmax .* max(0, 1 - (r_v ./ max(R_eff, eps)).^2);
v(r_v < eps) = Vmax;
end
