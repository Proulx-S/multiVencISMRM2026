% dev_signal_sanity.m -- verify I normalization for centerVox and pseudoVoxel
addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

p = runSim;
pS = p.pSim;
pS.voxGrid.fovFE = 3.1; pS.voxGrid.fovPE = 3;
pS.voxGrid.matFE = 3; pS.voxGrid.matPE = 5;
pS.nSpin = (2^5)^2;
pV = p.pVessel;
pV.ID = pS.voxGrid.fovPE .* 0.5;
pV.WT = pS.voxGrid.fovPE .* 1;
pV.S.lumen   = 2.6879;
pV.S.wall = 5.234;
pV.S.surround = 1.32;


res = runSim(pV, pS, [], [], false);
res.pSim.voxGrid
res.pSim.spinGrid
% imagesc(res.vMap)
% ax = gca; ax.DataAspectRatio = [1 1 1];
imagesc(res.pSim.spinGrid.coorPE,res.pSim.spinGrid.coorFE,res.magMap)
axis image
xline(res.pSim.voxGrid.coorPE(1:end-1)+res.pSim.voxGrid.dPE/2,'k')
yline(res.pSim.voxGrid.coorFE(1:end-1)+res.pSim.voxGrid.dFE/2,'k')
colorbar

[voxGridFE,voxGridPE] = ndgrid(res.pSim.voxGrid.coorFE, res.pSim.voxGrid.coorPE);
voxGridR = sqrt(voxGridFE.^2+voxGridPE.^2);
imagesc(res.pSim.voxGrid.coorPE,res.pSim.voxGrid.coorFE,sqrt(voxGridFE.^2+voxGridPE.^2))

voxIdx = getVoxIdx(res.pSim.voxGrid, res.pSim.spinGrid);
sum(res.magMap(voxIdx==max(unique(voxIdx))))
sum(res.magMap(voxIdx==0))

imagesc(voxIdx)
res.pVessel.S
[~,b] = min(voxGridR(:));
