% dev_matchedSimSummary.m -- matched simulation figure in same format as phantomSummary
% Run from: cd /scratch/bass/projects/multiVencISMRM2026/tmp && matlab -batch "run('dev_matchedSimSummary.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
saveThis = 1;

%% Phantom-matched vessel parameters
ID      =  6.35;    % mm
WT      =  2.38125; % mm
OD      =  ID + 2*WT;
vMean   =  4.234;   % cm/s  -- from radial profile fit

%% Fine M1-spaced venc list for smooth spiral
gamma_phys = 2.6752218708e8 / (2*pi); % Hz/T
M1_max     = vencToM1(2);             % M1 for venc=2 cm/s (strongest practical encoding)
nSteps     = 400;
M1_fine    = linspace(M1_max/nSteps, M1_max, nSteps); % evenly spaced in M1
venc_fine  = pi * 100 ./ (gamma_phys .* M1_fine);     % corresponding venc [cm/s]

%% Simulation parameters
pDef = runSim; % get defaults

% Vessel: matched to phantom
pVesselSim          = pDef.pVessel;
pVesselSim.ID       = ID;
pVesselSim.WT       = WT;
pVesselSim.vMean    = vMean;
pVesselSim.profile  = 'parabolic1';
pVesselSim.S.lumen  = []; % velocity-dependent inflow (matched physics)

% Spin grid: FOV large enough to include vessel + wall + some tissue
pSimSim          = pDef.pSim;
pSimSim.fovFE    = (OD + 3) * 1.0; % OD + 3mm margin  ≈ 16.9mm
pSimSim.fovPE    = (OD + 5) * 1.0; % slightly wider    ≈ 18.9mm
pSimSim.matFE    = 3;
pSimSim.matPE    = 3;
pSimSim.nSpin    = (2^8+1)^2; % 257^2 spins → fine spatial resolution

% MRI: PCmono with fine M1 sweep
pMriSim                  = pDef.pMri;
pMriSim.venc.method      = 'PCmono';
pMriSim.venc.vencList    = venc_fine(:);
pMriSim.venc.FVEres      = 0;   pMriSim.venc.FVEbw     = 0;
pMriSim.venc.FVEvel      = [];  pMriSim.venc.vencMin    = [];
pMriSim.venc.vencMax     = [];

fprintf('Running matched simulation (FOV=%.1fx%.1fmm, %d venc steps)...\n', ...
    pSimSim.fovFE, pSimSim.fovPE, nSteps);
res = runSim(pVesselSim, pSimSim, pMriSim, [], false); % light=false → keep magMap/vMap
fprintf('Done.\n');

%% Coordinate axes for spatial maps (spin grid positions)
FEax = res.pSim.gridFE(:,1); % FE axis [mm]
PEax = res.pSim.gridPE(1,:); % PE axis [mm]

%% Complex-plane spiral
Iref    = res.I(1,1,1,1,1,2);           % M1=0 reference signal (real, positive)
Ienc    = squeeze(res.I(1,1,1,1,:,1));  % velocity-encoded [nSteps x 1]
Mnorm   = abs(Iref);
Ienc_norm = Ienc / Mnorm;               % normalized complex trajectory

%% Build figure: same format as phantomSummary (tiledlayout 3x5, columnmajor)
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT = tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor');
ax = {};

% --- Tile 1: magnitude map ---
ax{end+1} = nexttile(hT);
imagesc(ax{end}, PEax, FEax, res.magMap, [0 max(res.magMap(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},'simulation ROI');

% --- Tile 2: velocity map ---
ax{end+1} = nexttile(hT);
vLim = max(abs(res.vMap(:)));
imagesc(ax{end}, PEax, FEax, res.vMap, [-vLim vLim]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = redblue; set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},sprintf('parabolic vMean=%.1f cm/s', vMean));

% --- Tile 3: lumen (blood-only) mask ---
ax{end+1} = nexttile(hT);
imagesc(ax{end}, PEax, FEax, single(res.pVessel.mask.lumen), [0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
theta_c = linspace(0,2*pi,360);
hold(ax{end},'on'); plot(ax{end}, ID/2*cos(theta_c), ID/2*sin(theta_c), 'm');
title(ax{end},'lumen mask');

% --- Tile 4: wall mask ---
ax{end+1} = nexttile(hT);
imagesc(ax{end}, PEax, FEax, single(res.pVessel.mask.wall), [0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
hold(ax{end},'on');
plot(ax{end}, ID/2*cos(theta_c), ID/2*sin(theta_c), 'm');
plot(ax{end}, OD/2*cos(theta_c), OD/2*sin(theta_c), 'm');
title(ax{end},'wall mask');

% --- Tile 5: tissue (surround) mask ---
ax{end+1} = nexttile(hT);
imagesc(ax{end}, PEax, FEax, single(res.pVessel.mask.surround), [0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
hold(ax{end},'on');
plot(ax{end}, OD/2*cos(theta_c), OD/2*sin(theta_c), 'm');
title(ax{end},'surround mask');

% --- Tiles [3x3]: complex-plane spiral (smooth colored line, fine M1 sweep) ---
ax{end+1} = nexttile(hT, [3 3]);

% Colored line using surface trick (colored by M1 / venc)
xSp = real(Ienc_norm);
ySp = imag(Ienc_norm);
cVal = M1_fine / M1_max; % normalized [0,1]: 0 = highest venc, 1 = lowest venc
hS = surface(ax{end}, [xSp(:) xSp(:)], [ySp(:) ySp(:)], ...
             zeros(nSteps,2), [cVal(:) cVal(:)], ...
             'EdgeColor','flat','FaceColor','none','LineWidth',1.5);
colormap(ax{end}, jet);
cb = colorbar(ax{end},'Location','eastoutside');
ylabel(cb, 'normalized M_1  (0=high venc, 1=low venc)');
hold(ax{end},'on');

% Reference circle (radius 1 in normalized units)
plot(ax{end}, cos(theta_c), sin(theta_c), 'w--', 'LineWidth', 0.8);
% Axes through origin
xline(ax{end},0,'w','LineWidth',0.5,'Alpha',0.5);
yline(ax{end},0,'w','LineWidth',0.5,'Alpha',0.5);

axis(ax{end},'image','tight');
xLim = xlim(ax{end}); if xLim(1)>0; xLim(1)=0; end; if xLim(2)<0; xLim(2)=0; end;
yLim = ylim(ax{end}); if yLim(1)>0; yLim(1)=0; end; if yLim(2)<0; yLim(2)=0; end;
dLim = max(diff(xLim),diff(yLim))*0.05;
set(ax{end},'XLim',xLim+[-dLim dLim],'YLim',yLim+[-dLim dLim]);
grid(ax{end},'on');
xlabel(ax{end},'real'); ylabel(ax{end},'imag');
ax{end}.Color = 'k';
ax{end}.GridColor = [0.5 0.5 0.5];
title(ax{end},'complex-domain signal evolution (simulation)');

drawnow;
if saveThis
    saveas(        f, fullfile(figDir, 'matchedSimSummary.fig'));
    exportgraphics(f, fullfile(figDir, 'matchedSimSummary.png'));
    exportgraphics(f, fullfile(figDir, 'matchedSimSummary.svg'));
    fprintf('Saved matchedSimSummary\n/local/users/Proulx-S/projects/multiVencISMRM2026/figures/matchedSimSummary.png\n');
end
disp('dev_matchedSimSummary done');
