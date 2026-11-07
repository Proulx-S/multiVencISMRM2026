function [Mz_n, pMri, Mo] = getMz_n(pMri, pRelax, n, Mo)
% getMz_n  Instantaneous longitudinal magnetization after exactly n RF pulses from fully-relaxed M0.
%
% Unlike getMz_ss (which averages over entry phases 1..n), getMz_n returns the
% magnetization of a spin that entered the slab just before the first pulse and
% has since received exactly n pulses. This is the correct model for isochromats
% at known Z positions within the slab (used by the 'cylinder3D' profile).
%
% NOTE (2026-11-06) -- version implemented here: CONTINUOUS COUNT.
%   Per isochromat Q1^n with a continuous (non-integer) n = d/(v*TR) passed in by simVesselSpins.
%   A spin only receives whole pulses, so a fractional n credits it with pulses it has not
%   received (Q1^x < Q1^floor(x)); the resulting inflow curve is ~12% low between integer n.
%   The exact version averages Q1^floor(x) over the isochromat's cell (piecewise linear in v;
%   column average f = [(1-Q1^K)/(1-Q1) + (U-K)*Q1^K]/U, K = floor(U)).
%   Refs: Gullberg, Simons & Wehrli 1988, Magn Reson Imaging 6(4):437-461,
%   doi:10.1016/0730-725X(88)90481-X (fractional-volume method); Gao, Holland & Gore 1988,
%   Med Phys 15(6):809-814, doi:10.1118/1.596197 (that method for spoiled FLASH);
%   Bianciardi et al. 2016, Phil Trans R Soc A 374:20150184, doi:10.1098/rsta.2015.0184.
%   Note: https://github.com/Proulx-S/microscopyValidation/blob/9057997d65a96cad63419fed24daa1ba885ed6c6/docs/inflowBetweenIntegerPulseCounts.md (private repo)
%   Path: no other getMz_n.m exists in util, pcMRAsim (dev-inflow), multiVencSim or vfMRItools,
%   so this copy is the one the doIt's path finds (nothing shadows it).
%
% INPUTS
%   pMri   - MRI acquisition struct. Required fields: TR, FA.
%   pRelax - Relaxation struct with field T1 [s]. Used to set (override) pMri.E1 and pMri.Q1.
%   n      - Number of RF pulses received (scalar or array, Inf → returns Mz_ss).
%   Mo     - (optional) Equilibrium magnetization. Default 1.
%
% OUTPUTS
%   Mz_n   - Longitudinal magnetization after n pulses. Same size as n.
%   pMri   - Same struct with E1 and Q1 updated.
%   Mo     - Equilibrium magnetization used.
%
% Formula: Mz_n = Mz_ss + (Mo - Mz_ss) * Q1^n
% where    Mz_ss = Mo*(1-E1)/(1-Q1),  Q1 = E1*cosd(FA),  E1 = exp(-TR/T1)
%
% Special cases:
%   n = 0    → Mz_n = Mo         (not yet pulsed, fully relaxed)
%   n = Inf  → Mz_n = Mz_ss      (fully saturated, stationary-spin steady state)

if exist('pRelax','var') && ~isempty(pRelax) && isfield(pRelax,'T1') && ~isempty(pRelax.T1)
    pMri.E1 = exp(-pMri.TR / pRelax.T1);
    pMri.Q1 = pMri.E1 * cosd(pMri.FA);
elseif ~isfield(pMri,'E1') || isempty(pMri.E1)
    error(['E1 not found in pMri' newline 'provide pRelax input']);
end
if ~exist('Mo','var') || isempty(Mo)
    Mo = 1;
end

Mz_ss = Mo * (1 - pMri.E1) / (1 - pMri.Q1);
Mz_n  = Mz_ss + (Mo - Mz_ss) .* pMri.Q1.^n;
end
