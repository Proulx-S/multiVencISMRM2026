% dev_slide08_matchedSim.m -- dev script for matched simulation overlay (slide 8)
% Run from: cd /scratch/bass/projects/multiVencISMRM2026 && matlab -batch "run('tmp/dev_slide08_matchedSim.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
saveThis = 1;

%% Load phantom data
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
data = conj(data); dataNoFlow = conj(dataNoFlow);

ID = 6.35; OD = 11.11; bestVenc = 9;

% Recompute masks (same as doIt.m)
M0 = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, size(M0,1)*FEspacing-FEspacing/2, size(M0,1));
PEpos = linspace(PEspacing/2, size(M0,2)*PEspacing-PEspacing/2, size(M0,2));
[FEgrid, PEgrid] = ndgrid(FEpos, PEpos);
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
maskWallLowMag = single(M0<0.44e-7);
theta_circ = linspace(0,2*pi,360);

%% Radial profile fit (recompute from block 3)
PD_tmp = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))), 3));
m1best = vencToM1(bestVenc);
velPD_tmp = phase2vel(PD_tmp, m1best);
if mean(velPD_tmp(maskBloodOnly)) < 0; velPD_tmp = -velPD_tmp; end
rFit3 = double(rGrid(maskBloodOnly)); vFit3 = double(velPD_tmp(maskBloodOnly));
parabola_fun3 = @(vMax, r) vMax .* (1 - (r./(ID/2)).^2);
opts3 = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000);
vMax_fit3 = lsqcurvefit(parabola_fun3, max(vFit3), rFit3(:), vFit3(:), 0, [], opts3);
vMean_fit3 = vMax_fit3 / 2;
fprintf('Fitted: vMax=%.3f cm/s, vMean=%.3f cm/s\n', vMax_fit3, vMean_fit3);

%% Run matched simulation at phantom venc levels
% Use PCmono matching the phantom acquisition
vencListPhantom = sort(unique(dataVenc(~isinf(dataVenc))),'descend'); % [112 56 ... 2] cm/s

pDefault = runSim; % get defaults

pVessel4 = pDefault.pVessel;
pVessel4.ID    = ID;          % 6.35mm, matched to phantom
pVessel4.WT    = 2.38125;     % wall thickness from tube spec
pVessel4.vMean = vMean_fit3;  % from radial profile fit
pVessel4.profile = 'parabolic1';
pVessel4.S.lumen = []; % velocity-dependent inflow (matches phantom conditions)

pSim4 = pDefault.pSim;
pSim4.fovFE = FEspacing; % match one phantom pixel
pSim4.fovPE = PEspacing;
pSim4.matFE = 3;
pSim4.matPE = 3;
pSim4.nSpin = (2^8+1)^2;

pMri4 = pDefault.pMri;
pMri4.venc.method   = 'PCmono';
pMri4.venc.vencList  = vencListPhantom(:);
pMri4.venc.FVEres = 0; pMri4.venc.FVEbw = 0;
pMri4.venc.FVEvel = []; pMri4.venc.vencMin = []; pMri4.venc.vencMax = [];

fprintf('Running matched simulation...\n');
res4 = runSim(pVessel4, pSim4, pMri4);
fprintf('Done.\n');

% Extract center voxel signal (lumen + surround)
% res4.I dims: [1 1 1 1 nVenc 2] = [FE PE SL t M1 M1ref]
Iref4 = res4.I(1,1,1,1,1,2);           % M1=0 reference (positive real after subtraction)
Ienc4 = squeeze(res4.I(1,1,1,1,:,1));  % velocity-encoded [nVenc x 1]
Isim4_all  = [Iref4; Ienc4(:)];
vsim4_all  = [inf; vencListPhantom(:)];
rsim4_all  = ones(length(Isim4_all),1);

%% Reproduce phantomSummary figure + simulated trajectory overlay
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT = tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};

% Mag
ax{end+1} = nexttile;
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
imagesc(ax{end},PEpos,FEpos,M,[0 max(M(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]); title(ax{end},'phantom ROI');

% Velocity map
ax{end+1} = nexttile;
PD_v   = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))),3));
velPD4 = phase2vel(PD_v, m1best);
imagesc(ax{end},PEpos,FEpos,velPD4,[-max(abs(velPD4(:))) max(abs(velPD4(:)))]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = redblue; set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},['venc=' num2str(bestVenc) 'cm/s']);

% Masks
ax{end+1} = nexttile; imagesc(ax{end},PEpos,FEpos,maskBloodOnly,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]); hold(ax{end},'on');
plot(ax{end},ID/2*cos(theta_circ), ID/2*sin(theta_circ), 'm'); title(ax{end},'blood-only mask');
ax{end+1} = nexttile; imagesc(ax{end},PEpos,FEpos,maskTissueOnly,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]); hold(ax{end},'on');
plot(ax{end},OD/2*cos(theta_circ), OD/2*sin(theta_circ), 'm'); title(ax{end},'tissue-only mask');
ax{end+1} = nexttile; imagesc(ax{end},PEpos,FEpos,maskWallLowMag,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]); hold(ax{end},'on');
plot(ax{end},ID/2*cos(theta_circ), ID/2*sin(theta_circ), 'm');
plot(ax{end},OD/2*cos(theta_circ), OD/2*sin(theta_circ), 'm'); title(ax{end},'low-mag wall mask');

% Complex domain: phantom data
ax{end+1} = nexttile([3 3]);
I = squeeze(mean(data,[1 2]));
Ivenc = squeeze(dataVenc);
Irun  = squeeze(dataRun);
[ax{end}, hP_phant] = plotMultiVenc(ax{end}, I, Ivenc, Irun, [], 'hot');
hold(ax{end},'on');

% Complex domain: simulated trajectory overlay (squares, black-outlined)
Mnorm4 = abs(Iref4); % reference magnitude
vencListSim_sorted = sort(unique(vsim4_all(~isinf(vsim4_all))),'descend');
M1sorted = vencToM1(vencListSim_sorted);
M1norm4  = M1sorted - min(M1sorted); M1norm4 = M1norm4 / max(M1norm4);
cMapSim4 = cool(length(vencListSim_sorted));
for kk = 1:length(vencListSim_sorted)
    v_ = vencListSim_sorted(kk);
    idx_ = vsim4_all == v_;
    IencNorm = real(Isim4_all(idx_)) / Mnorm4;
    IencImag = imag(Isim4_all(idx_)) / Mnorm4;
    plot(ax{end}, IencNorm, IencImag, 's', 'MarkerFaceColor', cMapSim4(kk,:), ...
         'MarkerEdgeColor', 'w', 'MarkerSize', 6);
end

title(ax{end},'complex-domain signal evolution (hot=phantom, cool=simulation)');

drawnow;
if saveThis
    saveas(        f, fullfile(figDir, 'phantomSummaryWithSim.fig'));
    exportgraphics(f, fullfile(figDir, 'phantomSummaryWithSim.png'));
    exportgraphics(f, fullfile(figDir, 'phantomSummaryWithSim.svg'));
    disp('Saved phantomSummaryWithSim');
end
disp('dev_slide08_matchedSim done');
