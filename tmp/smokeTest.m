addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/util'));
cd('/scratch/bass/projects/multiVencISMRM2026');

disp('=== getMz_n smoke test ===');
pMri.TR = 0.01267; pMri.FA = 50;
pRelax.T1 = 3.25;
[Mz0,  pMri] = getMz_n(pMri, pRelax, 0);
[Mzss, ~]    = getMz_n(pMri, [],     Inf);
% Analytical Mz_ss: Mo*(1-E1)/(1-Q1)
E1 = exp(-pMri.TR / pRelax.T1);
Q1 = E1 * cosd(pMri.FA);
Mzss_ref = (1 - E1) / (1 - Q1);
fprintf('Mz_n(n=0)=%.4f  Mz_n(n=Inf)=%.4f  Mz_ss_ref=%.4f\n', Mz0, Mzss, Mzss_ref);
assert(abs(Mz0 - 1) < 1e-10,          'n=0 should return Mo=1');
assert(abs(Mzss - Mzss_ref) < 1e-10,  'n=Inf should match Mz_ss');
disp('getMz_n: PASS');

disp('=== simCylinder3D smoke test ===');
p_def = runSim();
pV = p_def.pVessel; pS = p_def.pSim; pM = p_def.pMri;
pM.fieldStrength  = 3;    pM.species = 'phantom';
pM.sliceThickness = 2.2;  pM.TR = 0.01265; pM.TE = 0.0098; pM.FA = 50;
pM.venc.method = 'FVEmono'; pM.venc.FVEbw = 100;
pV.profile  = 'cylinder3D';
pV.ID       = 6.35;  pV.WT   = 2.38;  pV.Vmax = 8;
pV.A        = 0.9;   pV.n_hat = [0;0;1];
pV.posFE    = 0;     pV.posPE = 0;    pV.posSLC = 0;
pV.S.surround = 0.2;
pS.voxGrid.fovFE = 9; pS.voxGrid.fovPE = 9;
pS.voxGrid.matFE = 3; pS.voxGrid.matPE = 3;
pS.nSpin = 49; pS.monteCarloN = 0;
res = runSim(pV, pS, pM, false, false, true);
fprintf('magMap size: %dx%d, range: %.4f–%.4f\n', ...
    size(res.magMap,1), size(res.magMap,2), min(res.magMap(:)), max(res.magMap(:)));
assert(~any(isnan(res.magMap(:))), 'NaN in magMap');
disp('simCylinder3D: PASS');

disp('ALL SMOKE TESTS PASSED');
exit
