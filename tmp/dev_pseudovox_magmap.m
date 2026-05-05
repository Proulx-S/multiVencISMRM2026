% dev_pseudovox_magmap.m -- test pseudoVoxel grid mode, show magMap
% Run from: cd /scratch/bass/projects/multiVencISMRM2026 && matlab -batch "run('tmp/dev_pseudovox_magmap.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

figDir = '/scratch/bass/projects/multiVencISMRM2026/tmp';

%% Vessel parameters (matched to phantom tube)
pDef = runSim;

pVessel          = pDef.pVessel;
pVessel.ID       = 6.35;     % [mm] tube inner diameter
pVessel.WT       = 2.38;     % [mm] wall thickness
pVessel.vMean    = 4.23;     % [cm/s] from radial profile fit (block 3)
pVessel.profile  = 'parabolic1';
pVessel.S.lumen  = [];       % auto (velocity-dependent inflow)
pVessel.S.wall   = 0;

%% Simulation grid — pseudoVoxel mode, 3x3 voxels
pSim             = pDef.pSim;
pSim.voxGrid.fovFE = 19.05;    % [mm] 3 x 6.35mm (one vessel diameter per voxel)
pSim.voxGrid.fovPE = 19.05;
pSim.voxGrid.matFE = 3;
pSim.voxGrid.matPE = 3;
pSim.nSpin       = (2^8+1)^2;  % target total spins
pSim.gridMode    = 'pseudoVoxel';

%% MRI parameters
pMri             = pDef.pMri;
pMri.fieldStrength = 3;
pMri.species       = 'phantom';
pMri.venc.method   = 'PCmono';
pMri.venc.vencList = [40; 20];
pMri.venc.FVEres = 0; pMri.venc.FVEbw = 0;
pMri.venc.FVEvel = []; pMri.venc.vencMin = []; pMri.venc.vencMax = [];

%% Run
fprintf('Grid: %dx%d voxels, nSpin target=%d\n', pSim.voxGrid.matFE, pSim.voxGrid.matPE, pSim.nSpin);
res = runSim(pVessel, pSim, pMri, [], false);  % light=false → keep magMap/vMap
nSpinPerVox = res.pSim.nSpinPerVox;
fprintf('Grid actual: %dx%d spins total, %d per voxel, dFE=%.4fmm\n', ...
    res.pSim.spinGrid.matPE, res.pSim.spinGrid.matFE, nSpinPerVox, res.pSim.spinGrid.dFE);

%% Display magMap
magDisp   = double(res.magMap) .* nSpinPerVox;  % undo normalization
% dim1=FE (rows), dim2=PE (cols) — no transpose needed

FEax = res.pSim.spinGrid.coorFE;
PEax = res.pSim.spinGrid.coorPE;

f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 20 10]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

% magMap
ax1 = nexttile;
imagesc(ax1, PEax, FEax, magDisp);
axis(ax1,'image'); colormap(ax1,'gray'); colorbar(ax1);
xlabel(ax1,'PE [mm]'); ylabel(ax1,'FE [mm]');
title(ax1,sprintf('magMap (pseudoVoxel, %dx%d spins)', res.pSim.spinGrid.matPE, res.pSim.spinGrid.matFE));
% voxel grid lines
for k = res.pSim.voxGrid.coorFE(1:end-1) + res.pSim.voxGrid.dFE/2; xline(ax1, k, 'r', 'LineWidth', 1); end
for k = res.pSim.voxGrid.coorPE(1:end-1) + res.pSim.voxGrid.dPE/2; yline(ax1, k, 'r', 'LineWidth', 1); end
% vessel boundary
theta = linspace(0,2*pi,360);
hold(ax1,'on');
plot(ax1, pVessel.ID/2.*cos(theta), pVessel.ID/2.*sin(theta), 'c--', 'LineWidth',1.5);
plot(ax1, (pVessel.ID/2+pVessel.WT).*cos(theta), (pVessel.ID/2+pVessel.WT).*sin(theta), 'y--', 'LineWidth',1.5);

% vMap
ax2 = nexttile;
vMapDisp = double(res.vMap);
imagesc(ax2, PEax, FEax, vMapDisp);
axis(ax2,'image'); colorbar(ax2);
xlabel(ax2,'PE [mm]'); ylabel(ax2,'FE [mm]');
title(ax2,'vMap [cm/s]');
for k = res.pSim.voxGrid.coorFE(1:end-1) + res.pSim.voxGrid.dFE/2; xline(ax2, k, 'r', 'LineWidth', 1); end
for k = res.pSim.voxGrid.coorPE(1:end-1) + res.pSim.voxGrid.dPE/2; yline(ax2, k, 'r', 'LineWidth', 1); end
hold(ax2,'on');
plot(ax2, pVessel.ID/2.*cos(theta), pVessel.ID/2.*sin(theta), 'c--', 'LineWidth',1.5);
plot(ax2, (pVessel.ID/2+pVessel.WT).*cos(theta), (pVessel.ID/2+pVessel.WT).*sin(theta), 'y--', 'LineWidth',1.5);

drawnow;
outFile = fullfile(figDir, 'dev_pseudovox_magmap.png');
exportgraphics(f, outFile, 'Resolution', 150);
fprintf('Saved: %s\n', outFile);
