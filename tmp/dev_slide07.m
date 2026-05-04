% dev_slide07.m -- dev script for simulation summary (slide 7)
% Run from: cd /scratch/bass/projects/multiVencISMRM2026 && matlab -batch "run('tmp/dev_slide07.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

saveThis = 1;

%% Get default simulation parameters
p = runSim;

%% Vessel parameters: two flow profiles, both with FLAT magnitude profile
% Flat magnitude: S.lumen = constant Mxy at vMean (no velocity-dependent inflow saturation)
Mz_flat  = getMz_ss(p.pMri, p.pMri.relax.blood, p.pVessel.vMean);
Mxy_flat = getMxy_ss(Mz_flat, p.pMri, p.pMri.relax.blood);

% parabolic (laminar) flow: default PD=0
pVesselPara = p.pVessel;
pVesselPara.S.lumen = Mxy_flat; % flat: all spins same magnitude

% plug flow: PD = ID (full plug)
pVesselPlug = p.pVessel;
pVesselPlug.PD      = pVesselPlug.ID;
pVesselPlug.S.lumen = Mxy_flat; % flat: all spins same magnitude

%% MRI parameters: PCmono with venc list matching in vivo data range
p.pMri.venc.method   = 'PCmono';
vencListSim           = [40 20 13 10 8 7 6 5 4]'; % [cm/s] matching in vivo
p.pMri.venc.vencList  = vencListSim;
p.pMri.venc.FVEres    = 0;
p.pMri.venc.FVEbw     = 0;
p.pMri.venc.FVEvel    = [];
p.pMri.venc.vencMin   = [];
p.pMri.venc.vencMax   = [];

%% Simulation parameters: small nSpin for speed (enough for maps)
p.pSim.nSpin = (2^7+1)^2;

disp('Running parabolic (laminar) simulation...');
resPara = runSim(pVesselPara, p.pSim, p.pMri, [], false); % light=false → keep magMap/vMap
disp('Running plug flow simulation...');
resPlug = runSim(pVesselPlug, p.pSim, p.pMri, [], false);
disp('Simulations done.');

%% Extract complex-domain signal for plotMultiVenc
% res.I dims: [FE PE SL t M1 M1ref] = [3 3 1 1 nVenc 2]
% After PCmono ref subtraction: col 1 = velocity-encoded, col 2 = real-valued reference (M1=0)
% res.I dims: [1 1 1 1 nVenc 2] = [FE PE SL t M1 M1ref]
% (FE/PE=1 because only center voxel is kept)
% col 1 = velocity-encoded signal, col 2 = M1=0 reference (real after phase subtraction)

for flowIdx = 1:2
    if flowIdx == 1; res = resPara; else; res = resPlug; end

    Iref  = res.I(1,1,1,1,1,2); % M1=0 reference signal (real-valued after subtraction)
    Ienc  = squeeze(res.I(1,1,1,1,:,1)); % velocity-encoded [nVenc x 1]
    Iplot  = [Iref; Ienc(:)];
    vPlot  = [inf; vencListSim];
    runPlot = ones(length(Iplot),1);

    if flowIdx == 1
        IplotPara = Iplot; vPlotPara = vPlot; runPlotPara = runPlot;
    else
        IplotPlug = Iplot; vPlotPlug = vPlot; runPlotPlug = runPlot;
    end
end

%% Build figure: 2 rows (plug, laminar) × 3 cols (magMap | vMap | complex-plane spiral)
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 22]);
hT = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact'); ax = {};

flowNames = {'plug flow','laminar flow'};
resList   = {resPlug,   resPara};
IplotList = {IplotPlug, IplotPara};
vList     = {vPlotPlug, vPlotPara};
runList_  = {ones(length(IplotPlug),1), ones(length(IplotPara),1)};

for rowIdx = 1:2
    res_   = resList{rowIdx};

    % magMap
    ax{end+1} = nexttile(hT);
    imagesc(res_.magMap); axis image;
    ax{end}.Colormap = gray;
    ax{end}.CLim = [0 max(res_.magMap(:))*1.1];
    set(ax{end},'XTick',[],'YTick',[]);
    title(ax{end},{flowNames{rowIdx},'magnitude map'});

    % vMap (velocity map) with diverging colormap
    ax{end+1} = nexttile(hT);
    vLim = max(abs(res_.vMap(:)))*1.1;
    imagesc(res_.vMap,[-vLim vLim]); axis image;
    ax{end}.Colormap = redblue;
    set(ax{end},'XTick',[],'YTick',[]);
    ylabel(colorbar,'velocity (cm/s)','FontSize',8);
    title(ax{end},{flowNames{rowIdx},'velocity map'});

    % Complex-domain signal evolution
    ax{end+1} = nexttile(hT);
    plotMultiVenc(ax{end}, IplotList{rowIdx}, vList{rowIdx}, runList_{rowIdx}, 'tight', 'jet');
    title(ax{end},{flowNames{rowIdx},'complex-domain signal'});
end

drawnow;
if saveThis
    saveas(        f, fullfile(figDir, 'simSummary.fig'));
    exportgraphics(f, fullfile(figDir, 'simSummary.png'));
    exportgraphics(f, fullfile(figDir, 'simSummary.svg'));
    disp('Saved simSummary');
end
disp('dev_slide07 done');
