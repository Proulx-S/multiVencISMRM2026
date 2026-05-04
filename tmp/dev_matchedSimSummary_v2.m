% dev_matchedSimSummary_v2.m -- corrected matched simulation figure
% Fixes: (1) grid transposition, (2) magnitude calibration, (3) phantom-matched FOV

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
saveThis = 1;

%% Load phantom data for calibration
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile,'data','dataVenc','dataRun','PEspacing','FEspacing');
data = conj(data);

ID = 6.35; OD = 11.1125; bestVenc = 10; % updated

nFE = size(data,1); nPE = size(data,2); % 27, 18

% Phantom coordinates + masks
M0 = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, nFE*FEspacing-FEspacing/2, nFE);
PEpos = linspace(PEspacing/2, nPE*PEspacing-PEspacing/2, nPE);
[FEgrid,PEgrid] = ndgrid(FEpos,PEpos);
total = sum(M0(:));
com(1) = sum(FEgrid(:).*M0(:))/total; com(2) = sum(PEgrid(:).*M0(:))/total;
FEgrid = FEgrid-com(1); FEpos = FEpos-com(1);
PEgrid = PEgrid-com(2); PEpos = PEpos-com(2);
rGrid = sqrt(PEgrid.^2+FEgrid.^2);
dFE = abs(FEgrid); dPE = abs(PEgrid);
d_far  = sqrt((dFE+FEspacing/2).^2+(dPE+PEspacing/2).^2);
d_near = sqrt(max(0,dFE-FEspacing/2).^2+max(0,dPE-PEspacing/2).^2);
maskBloodOnly  = d_far  < ID/2;
maskTissueOnly = d_near > OD/2;

% Magnitude calibration from phantom
M_blood  = double(mean(M0(maskBloodOnly)));
M_tissue = double(mean(M0(maskTissueOnly)));
fprintf('Phantom: M_blood=%.4g, M_tissue=%.4g, ratio=%.2f\n', M_blood, M_tissue, M_blood/M_tissue);

%% Radial profile fit (from block 3)
PD_tmp = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))),3));
m1bv   = vencToM1(bestVenc);
velPD_tmp = phase2vel(PD_tmp, m1bv);
if mean(velPD_tmp(maskBloodOnly)) < 0; velPD_tmp = -velPD_tmp; end
rFit   = double(rGrid(maskBloodOnly)); vFit = double(velPD_tmp(maskBloodOnly));
pfit   = @(vMax,r) vMax.*(1-(r./(ID/2)).^2);
opts   = optimoptions('lsqcurvefit','Display','off');
vMax_fit  = lsqcurvefit(pfit, max(vFit), rFit(:), vFit(:), 0, [], opts);
vMean_fit = vMax_fit / 2;
fprintf('Fit: vMax=%.3f cm/s, vMean=%.3f cm/s\n', vMax_fit, vMean_fit);

%% Simulation setup
pDef = runSim;

% Vessel: matched to phantom
pVesselSim          = pDef.pVessel;
pVesselSim.ID       = ID;
pVesselSim.WT       = 2.38125;
pVesselSim.vMean    = vMean_fit;
pVesselSim.profile  = 'parabolic1';
pVesselSim.S.lumen  = []; % velocity-dependent inflow (auto-computed)
pVesselSim.S.wall   = 0;  % tube wall: near-zero signal

% FOV exactly matching phantom
pSimSim          = pDef.pSim;
pSimSim.fovFE    = nFE * FEspacing; % 27 * 0.5 = 13.5mm
pSimSim.fovPE    = nPE * PEspacing; % 18 * 0.8929 = 16.07mm
pSimSim.matFE    = 3;
pSimSim.matPE    = 3;
pSimSim.nSpin    = (2^8+1)^2; % high resolution

% Fine M1-spaced venc for spiral (400 steps, evenly spaced in M1)
gamma_phys   = 2.6752218708e8 / (2*pi); % Hz/T
M1_max       = vencToM1(2);             % M1 for venc=2 cm/s
nSteps       = 400;
M1_fine      = linspace(M1_max/nSteps, M1_max, nSteps);
venc_fine    = pi * 100 ./ (gamma_phys .* M1_fine);  % cm/s (decreasing)

pMriSim                  = pDef.pMri;
pMriSim.venc.method      = 'PCmono';
pMriSim.venc.vencList    = venc_fine(:);
pMriSim.venc.FVEres      = 0; pMriSim.venc.FVEbw = 0;
pMriSim.venc.FVEvel      = []; pMriSim.venc.vencMin = []; pMriSim.venc.vencMax = [];

fprintf('Running simulation (FOV %.1fx%.1f mm, %d M1 steps)...\n', ...
    pSimSim.fovFE, pSimSim.fovPE, nSteps);
res = runSim(pVesselSim, pSimSim, pMriSim, [], false);
fprintf('Done.\n');

%% Calibrate S.surround to match phantom magnitude ratio
% magMap from runSim is already S.lumen/nSpin per spin — work in these per-spin units
lumen_mask    = res.pVessel.mask.lumen;
wall_mask     = res.pVessel.mask.wall;
surround_mask = res.pVessel.mask.surround;

magMap = double(res.magMap); % per-spin values: S.lumen/nSpin in lumen, S.surround/nSpin elsewhere

S_lumen_per_spin_mean = mean(magMap(lumen_mask));
S_surround_auto_per_spin = mean(magMap(surround_mask));
ratio_sim    = S_lumen_per_spin_mean / S_surround_auto_per_spin;
ratio_target = M_blood / M_tissue;

fprintf('Auto ratio = %.2f, target = %.2f\n', ratio_sim, ratio_target);

% Correct surround per-spin value to match phantom ratio
S_surround_corrected_per_spin = S_lumen_per_spin_mean / ratio_target;
magMap(surround_mask) = S_surround_corrected_per_spin;
magMap(wall_mask)     = 0;

fprintf('Corrected ratio: %.2f (target: %.2f)\n', ...
    mean(magMap(lumen_mask)) / mean(magMap(surround_mask)), ratio_target);

%% Grid axes: rows=PE, cols=FE in runSim → need transpose for display
% After transpose: rows=FE (y-axis), cols=PE (x-axis) — matches phantom convention
FEax_sim = res.pSim.gridFE(1,:); % FE values along cols (1×nGrid)
PEax_sim = res.pSim.gridPE(:,1); % PE values along rows (nGrid×1)
magMap_T  = magMap';              % [nGrid_FE, nGrid_PE] after transpose
vMap_T    = double(res.vMap)';
mask_lumen_T    = lumen_mask';
mask_wall_T     = wall_mask';
mask_surround_T = surround_mask';

theta_c = linspace(0,2*pi,360);

%% Complex-plane spiral
Iref_sim = res.I(1,1,1,1,1,2);
Ienc_sim = squeeze(res.I(1,1,1,1,:,1));
Ienc_norm = Ienc_sim / abs(Iref_sim);

%% Build figure: same 3x5 tiledlayout as phantomSummary
fSim = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hTSim = tiledlayout(fSim,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor');
axSim = {};

% Tile 1: magnitude map (rows=FE=y, cols=PE=x)
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, magMap_T, [0 max(magMap_T(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end},'simulation ROI');

% Tile 2: velocity map
axSim{end+1} = nexttile(hTSim);
vLim_sim = max(abs(vMap_T(:)));
imagesc(axSim{end}, PEax_sim, FEax_sim, vMap_T, [-vLim_sim vLim_sim]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axSim{end}.Colormap = redblue; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end}, sprintf('parabolic vMean=%.1f cm/s', vMean_fit));

% Tile 3: lumen mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(mask_lumen_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on');
plot(axSim{end}, ID/2*cos(theta_c), ID/2*sin(theta_c), 'm');
title(axSim{end},'lumen mask');

% Tile 4: wall mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(mask_wall_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on');
plot(axSim{end}, ID/2*cos(theta_c), ID/2*sin(theta_c), 'm');
plot(axSim{end}, OD/2*cos(theta_c), OD/2*sin(theta_c), 'm');
title(axSim{end},'wall mask');

% Tile 5: surround mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(mask_surround_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on');
plot(axSim{end}, OD/2*cos(theta_c), OD/2*sin(theta_c), 'm');
title(axSim{end},'surround mask');

% Tiles [3x3]: complex-plane spiral (smooth colored line, fine M1 sweep)
axSim{end+1} = nexttile(hTSim, [3 3]);
xSp  = real(Ienc_norm);
ySp  = imag(Ienc_norm);
cVal = (1:nSteps)' / nSteps;
surface(axSim{end}, [xSp xSp], [ySp ySp], zeros(nSteps,2), [cVal cVal], ...
    'EdgeColor','flat','FaceColor','none','LineWidth',1.5);
colormap(axSim{end}, jet);
cb_s = colorbar(axSim{end},'Location','eastoutside');
ylabel(cb_s,'normalized M_1  (0=low → 1=high)');
hold(axSim{end},'on');
plot(axSim{end}, cos(theta_c), sin(theta_c), 'w--', 'LineWidth', 0.8);
xline(axSim{end},0,'w','LineWidth',0.5,'Alpha',0.4);
yline(axSim{end},0,'w','LineWidth',0.5,'Alpha',0.4);
axis(axSim{end},'image','tight');
xL = xlim(axSim{end}); if xL(1)>0; xL(1)=0; end; if xL(2)<0; xL(2)=0; end;
yL = ylim(axSim{end}); if yL(1)>0; yL(1)=0; end; if yL(2)<0; yL(2)=0; end;
dL = max(diff(xL),diff(yL))*0.05;
set(axSim{end},'XLim',xL+[-dL dL],'YLim',yL+[-dL dL]);
grid(axSim{end},'on'); axSim{end}.Color = 'k'; axSim{end}.GridColor = [0.5 0.5 0.5];
xlabel(axSim{end},'real'); ylabel(axSim{end},'imag');
title(axSim{end},'complex-domain signal evolution (simulation)');

drawnow;
if saveThis
    saveas(        fSim, fullfile(figDir,'matchedSimSummary.fig'));
    exportgraphics(fSim, fullfile(figDir,'matchedSimSummary.png'));
    exportgraphics(fSim, fullfile(figDir,'matchedSimSummary.svg'));
    fprintf('Saved\n/local/users/Proulx-S/projects/multiVencISMRM2026/figures/matchedSimSummary.png\n');
end
disp('dev_matchedSimSummary_v2 done');
