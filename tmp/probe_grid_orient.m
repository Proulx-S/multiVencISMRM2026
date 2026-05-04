% Probe grid orientation and magnitude calibration
addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

% Phantom parameters
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile,'data','dataVenc','dataRun','dataNoFlow','PEspacing','FEspacing');
data = conj(data);
ID=6.35; OD=11.1125; bestVenc=10; % new bestVenc

nFE = size(data,1); nPE = size(data,2);
M0 = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));

% Phantom coordinate system
FEpos = linspace(FEspacing/2, nFE*FEspacing-FEspacing/2, nFE);
PEpos = linspace(PEspacing/2, nPE*PEspacing-PEspacing/2, nPE);
[FEgrid,PEgrid]=ndgrid(FEpos,PEpos);
total=sum(M0(:)); com(1)=sum(FEgrid(:).*M0(:))/total; com(2)=sum(PEgrid(:).*M0(:))/total;
FEgrid=FEgrid-com(1); FEpos=FEpos-com(1); PEgrid=PEgrid-com(2); PEpos=PEpos-com(2);
rGrid=sqrt(PEgrid.^2+FEgrid.^2);
dFE=abs(FEgrid); dPE=abs(PEgrid);
d_far=sqrt((dFE+FEspacing/2).^2+(dPE+PEspacing/2).^2);
d_near=sqrt(max(0,dFE-FEspacing/2).^2+max(0,dPE-PEspacing/2).^2);
maskBloodOnly=d_far<ID/2; maskTissueOnly=d_near>OD/2;

% Magnitude calibration targets
M_blood  = double(mean(M0(maskBloodOnly)));
M_tissue = double(mean(M0(maskTissueOnly)));
maskWallLowMag = single(M0 < 0.44e-7);
M_wall_region  = mean(M0(maskWallLowMag==1));

fprintf('Phantom FOV: FE=%.2fmm (%d pix x %.4fmm), PE=%.2fmm (%d pix x %.4fmm)\n', ...
    nFE*FEspacing, nFE, FEspacing, nPE*PEspacing, nPE, PEspacing);
fprintf('M_blood (blood-only pixels): %.4g\n', M_blood);
fprintf('M_tissue (tissue-only pixels): %.4g\n', M_tissue);
fprintf('M_wall (low-mag wall pixels): %.4g\n', M_wall_region);
fprintf('Ratio M_blood/M_tissue: %.2f\n', M_blood/M_tissue);

% Run a test simulation and check grid orientation
pDef = runSim;
pDef.pSim.fovFE = nFE*FEspacing;  % match phantom exactly
pDef.pSim.fovPE = nPE*PEspacing;
pDef.pSim.matFE = 3; pDef.pSim.matPE = 3;
pDef.pSim.nSpin = (2^6+1)^2; % small for speed
pDef.pVessel.ID = ID; pDef.pVessel.WT = 2.38125; pDef.pVessel.vMean = 4.234;
pDef.pMri.venc.method='PCmono'; pDef.pMri.venc.vencList=[10];
pDef.pMri.venc.FVEres=0; pDef.pMri.venc.FVEbw=0; pDef.pMri.venc.FVEvel=[];
pDef.pMri.venc.vencMin=[]; pDef.pMri.venc.vencMax=[];
res = runSim(pDef.pVessel, pDef.pSim, pDef.pMri, [], false);

gFE = res.pSim.gridFE; gPE = res.pSim.gridPE;
fprintf('\nGrid size: %dx%d\n', size(gFE,1), size(gFE,2));
fprintf('gridFE(:,1) range: [%.2f, %.2f] mm (should be FE = [%.2f, %.2f])\n', ...
    min(gFE(:,1)), max(gFE(:,1)), min(FEpos), max(FEpos));
fprintf('gridFE(1,:) range: [%.2f, %.2f] mm\n', min(gFE(1,:)), max(gFE(1,:)));
fprintf('gridPE(:,1) range: [%.2f, %.2f] mm\n', min(gPE(:,1)), max(gPE(:,1)));
fprintf('gridPE(1,:) range: [%.2f, %.2f] mm (should be PE = [%.2f, %.2f])\n', ...
    min(gPE(1,:)), max(gPE(1,:)), min(PEpos), max(PEpos));

% Check magMap: blood region should be in center, verify
fprintf('\nmagMap size: %dx%d\n', size(res.magMap,1), size(res.magMap,2));
fprintf('magMap range: [%.4g, %.4g]\n', min(res.magMap(:)), max(res.magMap(:)));
fprintf('S.lumen (auto-computed): %.4g\n', max(res.pVessel.S.lumen(:)));
fprintf('S.surround: %.4g\n', max(res.pVessel.S.surround(:)));
fprintf('S.wall: %g\n', res.pVessel.S.wall);
fprintf('Ratio S.lumen/S.surround: %.2f\n', mean(res.pVessel.S.lumen(:))/res.pVessel.S.surround);
fprintf('Target phantom ratio: %.2f\n', M_blood/M_tissue);
