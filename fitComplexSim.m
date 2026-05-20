function [fitResult, fitInfo] = fitComplexSim(s_meas, finiteVencs, pVessel_0, pSim_0, pMri, theta0, nSpinsSim)
% Fit vessel parameters to multi-VENC complex signal using runSim as forward model.
%
%   s_meas      — (nFE × nPE × K) complex, phase-corrected (ALL ROI pixels)
%   finiteVencs — (1×K) cm/s
%   pVessel_0   — base vessel struct (WT, posFE, posPE etc.; geometry overridden per iteration)
%   pSim_0      — base sim struct (voxGrid already set to ROI dimensions)
%   pMri        — resolved MRI params (relax.blood, relax.GM, TR, FA, etc.)
%   theta0      — [Vmax, R, e1, e2, A_n]
%                   e1 = (AR-1)*cos(2*alpha), e2 = (AR-1)*sin(2*alpha)
%                   A_n: normalised amplitude (predicted/measured scale factor, ~1 at init)
%   nSpinsSim   — target total spins (default 400 for speed)
%
%   Forward model: runSim with gridMode='allVoxels' returns per-pixel complex signal
%   at each M1 (VENC). S.lumen per spin = getMxy_ss(v_spin, blood T1) — physics-based.
%   Partial volume and within-voxel phase dispersion are handled by the spin summation.

if ~exist('nSpinsSim','var') || isempty(nSpinsSim); nSpinsSim = 400; end

[nFE, nPE, K] = size(s_meas);
s_scale = mean(abs(s_meas(:)));
s_n     = s_meas / s_scale;         % normalised measured signal (mean |s_n| ≈ 1)
s_std   = max(std(abs(s_n(:))), eps);

% --- Pre-compute spin grid (one-time cost, outside optimizer) ---
pSim_fit = pSim_0;
pSim_fit.nSpin       = nSpinsSim;
pSim_fit.gridMode    = 'allVoxels';
pSim_fit.monteCarloN = 0;

pV_setup = pVessel_0;
[AR_init, alpha_init] = e_to_ar_alpha(theta0(3), theta0(4));
pV_setup.ID = theta0(2) * 2;  pV_setup.AR = AR_init;  pV_setup.alpha = alpha_init;
pV_setup.vMean = theta0(1)/2;  pV_setup.vMax = [];  pV_setup.profile = 'parabolic1';
pV_setup.S.lumen = [];  pV_setup.S.surround = [];
% earlyStop=true: runs simVesselSpins only (sets up spinGrid, nSpinPerVox)
p_setup  = runSim(pV_setup, pSim_fit, pMri, false, false, true);
pSim_fit = p_setup.pSim;
% pMri is already resolved; use its m1List for VENC mapping
m1List  = double(pMri.venc.m1List(:,1));
m1_meas = arrayfun(@vencToM1, finiteVencs);
venc_idx = zeros(1,K);
for k = 1:K; [~, venc_idx(k)] = min(abs(m1List - m1_meas(k))); end

% Spin grid coordinates (FE, PE relative to grid centre; posFE=posPE=0 fixed)
[spinFEg, spinPEg] = ndgrid(pSim_fit.spinGrid.coorFE, pSim_fit.spinGrid.coorPE);

% Pre-compute tissue steady-state Mxy (constant across iterations)
Mxy_tissue = double(getMxy_ss(getMz_ss(pMri, pMri.relax.GM, 0), pMri, pMri.relax.GM));

% --- Bounds: internal theta = [Vmax, R, e1, e2, A_n] ---
% (e1,e2) = (AR-1)*[cos(2*alpha), sin(2*alpha)] — Cartesian ellipse coords
R_max  = min(pSim_fit.voxGrid.fovFE, pSim_fit.voxGrid.fovPE) / 2;
lb     = [0,   1e-6, -1, -1,  0];
ub     = [inf, R_max,  1,  1, inf];

% Initialise A_n
Mxy_ref  = double(getMxy_ss(getMz_ss(pMri,pMri.relax.blood,theta0(1)/2),pMri,pMri.relax.blood));
A_n_init = 1 / max(Mxy_ref, eps);
if numel(theta0) < 5 || isempty(theta0(5))
    theta0(5) = A_n_init;
end
% theta0 is already in (e1, e2) space
theta0 = min(max(theta0(:)', lb), ub);

% --- Optimizer ---
opts = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',2000);
theta = lsqnonlin( ...
    @(th) residuals_sim(th, spinFEg, spinPEg, pVessel_0, pSim_fit, pMri, ...
                        s_n, venc_idx, s_std, Mxy_tissue, nFE, nPE), ...
    theta0, lb, ub, opts);

% --- Output ---
[AR_f, alpha_f] = e_to_ar_alpha(theta(3), theta(4));
fitResult.Vmax     = theta(1);
fitResult.R        = theta(2);
fitResult.e1       = theta(3);
fitResult.e2       = theta(4);
fitResult.AR       = AR_f;     % derived
fitResult.alpha    = alpha_f;  % derived
fitResult.A_n      = theta(5);
fitResult.A        = theta(5) * s_scale;
fitResult.FEoffset = 0;
fitResult.PEoffset = 0;

[AR_0, alpha_0] = e_to_ar_alpha(theta0(3), theta0(4));
fitInfo.names  = {'Vmax','R','e1','e2','A_n','FEoffset','PEoffset'};
fitInfo.units  = {'cm/s','mm','-','-','-','mm','mm'};
fitInfo.fixed  = [false, false, false, false, false, true, true];
fitInfo.theta0 = [theta0(1:2), theta0(3:4), theta0(5), 0, 0];
fitInfo.lb     = [lb, 0, 0];
fitInfo.ub     = [ub, 0, 0];
fitInfo.theta  = [theta(1:2), theta(3:4), theta(5), 0, 0];
fitInfo.derived.names  = {'AR', 'alpha'};
fitInfo.derived.units  = {'-', 'rad'};
fitInfo.derived.theta0 = [AR_0, alpha_0];
fitInfo.derived.theta  = [AR_f, alpha_f];
fitInfo.nSpinsSim = nSpinsSim;
end


function res = residuals_sim(theta, spinFEg, spinPEg, pVessel_0, pSim_fit, pMri, ...
                              s_n, venc_idx, s_std, Mxy_tissue, nFE, nPE)
% One runSim call per iteration. Returns [real; imag] flattened residuals.
% theta = [Vmax, R, e1, e2, A_n]  (e1,e2 = Cartesian ellipse coords)
Vmax = theta(1); R = theta(2); A_n = theta(5);
[AR, alpha] = e_to_ar_alpha(theta(3), theta(4));

% Velocity at each spin (elliptical parabolic profile, centre at origin)
v_spins = velfunc_ellipse(spinFEg, spinPEg, Vmax, R, AR, alpha);

% Ellipse mask (analytical, same formula as simVesselSpins.m)
dFE = spinFEg;  dPE = spinPEg;
if abs(AR - 1) < 1e-6
    mask_lumen = sqrt(dFE.^2 + dPE.^2) <= R;
else
    u  =  dPE.*cos(alpha) + dFE.*sin(alpha);
    w  = -dPE.*sin(alpha) + dFE.*cos(alpha);
    mask_lumen = sqrt(u.^2 + (AR.*w).^2) <= R;
end

% Per-spin blood Mxy (physics, in normalised units)
Mxy_blood = double(getMxy_ss(getMz_ss(pMri,pMri.relax.blood,v_spins(mask_lumen)), ...
                              pMri, pMri.relax.blood));

% Build pVessel for this iteration
pV = pVessel_0;
pV.ID          = R * 2;
pV.AR          = AR;
pV.alpha       = alpha;
pV.vMean       = [];    % velocity supplied via profile
pV.vMax        = [];
pV.profile     = v_spins;               % numeric per-spin velocity [cm/s]
pV.S.lumen     = single(A_n .* Mxy_blood);   % per lumen-spin signal (normalised units)
pV.S.surround  = single(A_n * Mxy_tissue);   % scalar tissue signal
pV.S.wall      = 0;

% Run simulation (light=true: discard magMap/vMap/spinMap, keep res.I only)
p_sim = runSim(pV, pSim_fit, pMri, false, true, false);

% Extract predicted signal at measured VENCs: res.I is (nFE × nPE × 1 × 1 × nM1 × nM1ref)
K        = numel(venc_idx);
s_pred_n = zeros(nFE, nPE, K);
for k = 1:K
    s_pred_n(:,:,k) = double(p_sim.I(:,:,1,1,venc_idx(k),1));
end

diff = s_n - s_pred_n;
res  = [real(diff(:)); imag(diff(:))] / s_std;
end


function [e1, e2] = ar_alpha_to_e(AR, alpha)
d  = AR - 1;
e1 = d * cos(2*alpha);
e2 = d * sin(2*alpha);
end

function [AR, alpha] = e_to_ar_alpha(e1, e2)
d     = sqrt(e1^2 + e2^2);
AR    = 1 + d;
alpha = atan2(e2, e1) / 2;
end

function v = velfunc_ellipse(FEg, PEg, Vmax, R, AR, alpha)
% Elliptical parabolic velocity; peak at (FE=0, PE=0).
r_v   = sqrt(PEg.^2 + FEg.^2);
uPE   = PEg ./ max(r_v, eps);
uFE   = FEg ./ max(r_v, eps);
A_    = uPE.*cos(alpha) + uFE.*sin(alpha);
B_    = -uPE.*sin(alpha) + uFE.*cos(alpha);
R_eff = R ./ sqrt(A_.^2 + AR.^2 .* B_.^2);
v     = Vmax .* max(0, 1 - (r_v ./ max(R_eff, eps)).^2);
v(r_v < eps) = Vmax;
end
