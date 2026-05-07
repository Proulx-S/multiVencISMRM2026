function magBoundaryFit = fitMagBoundaryProfile(r, p, offset, m, rWO0, dr0, pLow0, pHigh0)
% Fit piecewise linear model to the amplitude radial profile of non-blood voxels.
%
%   r      — Nx1 radial distance from grid centre [mm]; rGrid(mask)
%   p      — Nx1 polar angle [rad]; pGrid(mask)
%             p = -atan2(FEgrid, PEgrid); p=0 → +PE axis (right in imagesc display)
%   offset — 1×2 [FEoffset, PEoffset] [mm]; vessel centre offset from grid centre
%   m      — Nx1 magnitudes of non-blood pixels (wall + tissue)
%   rWO0   — scalar; starting value for outer wall limit rWO [mm]
%   dr0    — scalar; starting value for transition width dr [mm] (e.g. mean voxel size)
%   pLow0  — scalar; starting value for wall plateau (e.g. mean wall-only magnitude)
%   pHigh0 — scalar; starting value for tissue plateau (e.g. mean tissue-only magnitude)
%
%   Model of offset-corrected radius r_off:
%     m(r_off) = pLow                                            r_off < rWO - dr/2
%     m(r_off) = pLow + (pHigh-pLow)*(r_off-(rWO-dr/2))/dr    rWO-dr/2 ≤ r_off ≤ rWO+dr/2
%     m(r_off) = pHigh                                           r_off > rWO + dr/2
%   Parameters: pLow (wall plateau), pHigh (tissue plateau), rWO (transition centre),
%               dr (transition width)
%
%   magBoundaryFit(r_off)   — cfit; input is offset-corrected radius [mm]
%   magBoundaryFit.rWO      — fitted outer wall limit [mm]
%   magBoundaryFit.dr       — fitted transition width [mm]
%   magBoundaryFit.pLow     — fitted wall signal level
%   magBoundaryFit.pHigh    — fitted tissue signal level

r = double(r(:));  p = double(p(:));  m = double(m(:));
FEoffset = double(offset(1));  PEoffset = double(offset(2));

% Normalise magnitudes so pHigh_init = 1 → coefficients O(1) regardless of signal scale
mScale = double(pHigh0);
m      = m      / mScale;
pLow0n = double(pLow0)  / mScale;
pHigh0n = 1;

% Offset-corrected radius: r_off = |pixel - vessel centre|
% with p = -atan2(FE, PE): r*cos(p) = PE, -r*sin(p) = FE
r_off = sqrt((r.*cos(p) - PEoffset).^2 + (-r.*sin(p) - FEoffset).^2);

% Reparameterise as [pLow, delta, rWO, dr] with delta = pHigh - pLow.
% Box constraints lb >= 0 then enforce pLow > 0, delta > 0 (i.e. pHigh > pLow), rWO > 0, dr > 0.
delta0n = pHigh0n - pLow0n;  % guaranteed > 0 since pHigh0 > pLow0 (tissue > wall)
wallModel = @(r_in, th) ...
    th(1) .* (r_in < th(3) - th(4)/2) + ...
    (th(1) + th(2) .* (r_in - (th(3)-th(4)/2)) ./ th(4)) ...
        .* (r_in >= th(3)-th(4)/2 & r_in <= th(3)+th(4)/2) + ...
    (th(1) + th(2)) .* (r_in > th(3) + th(4)/2);

th0  = double([pLow0n, delta0n, rWO0, dr0]);
lb   = [0,   0,   0,       0.01];
ub   = [inf, inf, 2*rWO0,  10  ];
opts = optimoptions('lsqnonlin', 'Display', 'off');
th   = lsqnonlin(@(th_) m - wallModel(r_off, th_), th0, lb, ub, opts);

% Convert [pLow, delta, rWO, dr] → [pLow, pHigh, rWO, dr] and un-normalise plateaus
th = [th(1)*mScale, (th(1)+th(2))*mScale, th(3), th(4)];

% Build cfit output — independent variable is offset-corrected radius
ft_boundary = fittype( ...
    @(pLow, pHigh, rWO, dr, r_in) ...
        pLow  .* (r_in < rWO - dr/2) + ...
        (pLow + (pHigh-pLow) .* (r_in-(rWO-dr/2)) ./ dr) ...
            .* (r_in >= rWO-dr/2 & r_in <= rWO+dr/2) + ...
        pHigh .* (r_in > rWO + dr/2), ...
    'independent', 'r_in', ...
    'coefficients', {'pLow', 'pHigh', 'rWO', 'dr'});

magBoundaryFit = cfit(ft_boundary, th(1), th(2), th(3), th(4));
end
