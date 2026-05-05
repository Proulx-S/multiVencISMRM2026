% dev_pseudovox_aniso.m -- test pseudoVoxel with anisotropic voxels and even matFE/matPE
% Run from: cd /scratch/bass/projects/multiVencISMRM2026 && matlab -batch "run('tmp/dev_pseudovox_aniso.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

figDir = '/scratch/bass/projects/multiVencISMRM2026/tmp';

pDef = runSim;

pVessel         = pDef.pVessel;
pVessel.ID      = 6.35;
pVessel.WT      = 2.38;
pVessel.vMean   = 4.23;
pVessel.profile = 'parabolic1';
pVessel.S.lumen = [];
pVessel.S.wall  = 0;

pSim           = pDef.pSim;
pSim.voxGrid.matFE = 4;
pSim.voxGrid.matPE = 5;
pSim.voxGrid.fovFE = pSim.voxGrid.matFE * 6.35;   % voxSzFE = 6.35 mm
pSim.voxGrid.fovPE = pSim.voxGrid.matPE * 4.00;   % voxSzPE = 4.00 mm  (anisotropic)
pSim.nSpin     = (2^8+1)^2;
pSim.gridMode  = 'pseudoVoxel';

pMri                  = pDef.pMri;
pMri.fieldStrength    = 3;
pMri.species          = 'phantom';
pMri.venc.method      = 'PCmono';
pMri.venc.vencList    = [40; 20];
pMri.venc.FVEres = 0; pMri.venc.FVEbw = 0;
pMri.venc.FVEvel = []; pMri.venc.vencMin = []; pMri.venc.vencMax = [];

fprintf('voxSzFE=%.2fmm  voxSzPE=%.2fmm  matFE=%d  matPE=%d\n', ...
    pSim.voxGrid.fovFE/pSim.voxGrid.matFE, pSim.voxGrid.fovPE/pSim.voxGrid.matPE, pSim.voxGrid.matFE, pSim.voxGrid.matPE);

res = runSim(pVessel, pSim, pMri, [], false);

nSpinPerVox = res.pSim.nSpinPerVox;
fprintf('Grid: %dx%d total spins | %d per voxel | dFE=%.4fmm dPE=%.4fmm\n', ...
    res.pSim.spinGrid.matPE, res.pSim.spinGrid.matFE, nSpinPerVox, res.pSim.spinGrid.dFE, res.pSim.spinGrid.dPE);

% Check boundary condition: voxel boundaries should fall exactly between spins
FEax = res.pSim.spinGrid.coorFE;  % 1 x nTotalFE
PEax = res.pSim.spinGrid.coorPE;  % 1 x nTotalPE

voxIdx = getVoxIdx(res.pSim.voxGrid, res.pSim.spinGrid);

% PE boundaries: vary along cols (dim 2) with ndgrid, fix middle FE row
gvPE = voxIdx(ceil(end/2), :);
PEbdry_actual = [];
for k = 1:length(gvPE)-1
    if gvPE(k) ~= gvPE(k+1)
        PEbdry_actual(end+1) = (PEax(k)+PEax(k+1))/2; %#ok<AGROW>
    end
end
voxBdryPE_theory = res.pSim.voxGrid.coorPE(1:end-1) + res.pSim.voxGrid.dPE/2;
fprintf('Theoretical PE boundaries: '); fprintf('%.4f ', voxBdryPE_theory); fprintf('mm\n');
fprintf('Actual    PE boundaries:   '); fprintf('%.4f ', PEbdry_actual);    fprintf('mm\n');
fprintf('Max PE boundary error: %.2e mm\n', max(abs(PEbdry_actual(:) - voxBdryPE_theory(:))));

% FE boundaries: vary along rows (dim 1) with ndgrid, fix middle PE col
gvFE = voxIdx(:, ceil(end/2));
FEbdry_actual = [];
for k = 1:length(gvFE)-1
    if gvFE(k) ~= gvFE(k+1)
        FEbdry_actual(end+1) = (FEax(k)+FEax(k+1))/2; %#ok<AGROW>
    end
end
voxBdryFE_theory = res.pSim.voxGrid.coorFE(1:end-1) + res.pSim.voxGrid.dFE/2;
fprintf('Theoretical FE boundaries: '); fprintf('%.4f ', voxBdryFE_theory); fprintf('mm\n');
fprintf('Actual    FE boundaries:   '); fprintf('%.4f ', FEbdry_actual);    fprintf('mm\n');
fprintf('Max FE boundary error: %.2e mm\n', max(abs(FEbdry_actual(:) - voxBdryFE_theory(:))));

%% Display
magDisp   = double(res.magMap) .* nSpinPerVox;
% dim1=FE (rows), dim2=PE (cols) — no transpose needed
vMapDisp = double(res.vMap);

f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 12]);
tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

theta = linspace(0,2*pi,360);

for iAx = 1:2
    ax = nexttile;
    if iAx==1
        imagesc(ax, PEax, FEax, magDisp); colormap(ax,'gray');
        title(ax, sprintf('magMap | %dx%d spins | dFE=%.3f dPE=%.3fmm', ...
            res.pSim.spinGrid.matPE, res.pSim.spinGrid.matFE, res.pSim.spinGrid.dFE, res.pSim.spinGrid.dPE));
    else
        imagesc(ax, PEax, FEax, vMapDisp);
        title(ax, 'vMap [cm/s]');
    end
    axis(ax,'image'); colorbar(ax);
    xlabel(ax,'PE [mm]'); ylabel(ax,'FE [mm]');
    for k = res.pSim.voxGrid.coorFE(1:end-1) + res.pSim.voxGrid.dFE/2; yline(ax, k, 'r', 'LineWidth', 1); end
    for k = res.pSim.voxGrid.coorPE(1:end-1) + res.pSim.voxGrid.dPE/2; xline(ax, k, 'r', 'LineWidth', 1); end
    hold(ax,'on');
    plot(ax, pVessel.ID/2.*cos(theta),              pVessel.ID/2.*sin(theta),              'c--','LineWidth',1.5);
    plot(ax,(pVessel.ID/2+pVessel.WT).*cos(theta),(pVessel.ID/2+pVessel.WT).*sin(theta),'y--','LineWidth',1.5);
end

drawnow;
outFile = fullfile(figDir,'dev_pseudovox_aniso.png');
exportgraphics(f, outFile, 'Resolution',150);
fprintf('Saved: %s\n', outFile);
