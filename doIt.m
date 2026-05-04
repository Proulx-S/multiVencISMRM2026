clear all; close all; clc;
% figure('MenuBar','none','ToolBar','none');


projectName = 'multiVencISMRM2026';
%%%%%%%%%%%%%%%%%%%%%
%% Set up environment
%%%%%%%%%%%%%%%%%%%%%

% Detect computing environment
os   = char(java.lang.System.getProperty('os.name'));
host = char(java.net.InetAddress.getLocalHost.getHostName);
user = char(java.lang.System.getProperty('user.name'));

% Setup folders
if strcmp(os,'Linux') && strcmp(host,'takoyaki') && strcmp(user,'sebp')
    envId = 1;
    storageDrive  = '/local/users/Proulx-S/';
    scratchDrive  = '/scratch/bass/';
    databaseDrive        = fullfile(storageDrive, 'db'       );
    databasePhantomDrive = fullfile(storageDrive, 'dbPhantom');
    projectCode    = fullfile(scratchDrive, 'projects', projectName       ); if ~exist(projectCode           ,'dir'); mkdir(projectCode           ); end
    projectStorage = fullfile(storageDrive, 'projects', projectName       ); if ~exist(projectStorage        ,'dir'); mkdir(projectStorage        ); end
    projectScratch = fullfile(scratchDrive, 'projects', projectName, 'tmp'); if ~exist(projectScratch        ,'dir'); mkdir(projectScratch        ); end
    projectDataBase        = fullfile(databaseDrive                       ); if ~exist(projectDataBase       ,'dir'); mkdir(projectDataBase       ); end
    projectDataBasePhantom = fullfile(databasePhantomDrive                ); if ~exist(projectDataBasePhantom,'dir'); mkdir(projectDataBasePhantom); end
    toolDir        = fullfile(scratchDrive, 'tools'                       ); if ~exist(toolDir               ,'dir'); mkdir(toolDir               ); end
else
    envId = 2;

    mountPoint = '/Users/sebastienproulx/remote/takoyakiLocal';
    [status, ~] = system(['mount | grep ' mountPoint]);
    if status ~= 0
      system(['sshfs takoyaki:/local/users/Proulx-S ' mountPoint ' -o follow_symlinks,reconnect,allow_other']);
    end
    
    storageDrive   = '/Users/sebastienproulx/bass';
    scratchDrive   = '/Users/sebastienproulx/bass';
    databaseDrive = fullfile(mountPoint, 'db');
    databasePhantomDrive = fullfile(mountPoint, 'dbPhantom');
    projectCode     = fullfile(scratchDrive, 'projects', projectName);        if ~exist(projectCode    ,'dir'); mkdir(projectCode    ); end
    projectStorage  = fullfile(storageDrive, 'projects', projectName);        if ~exist(projectStorage ,'dir'); mkdir(projectStorage ); end
    projectScratch  = fullfile(scratchDrive, 'projects', projectName, 'tmp'); if ~exist(projectScratch ,'dir'); mkdir(projectScratch ); end
    projectDataBase        = fullfile(databaseDrive                        ); if ~exist(projectDataBase,'dir'); mkdir(projectDataBase); end
    projectDataBasePhantom = fullfile(databasePhantomDrive                 ); if ~exist(projectDataBasePhantom,'dir'); mkdir(projectDataBasePhantom); end
    toolDir         = fullfile(scratchDrive, 'tools'                       ); if ~exist(toolDir        ,'dir'); mkdir(toolDir        ); end
end

% Load dependencies and set paths
%%% initial cloning of matlab util to get gitClone.m
tool = 'util'; toolURL = 'https://github.com/Proulx-S/util.git';
if ~exist(fullfile(toolDir, tool), 'dir'); system(['git clone ' toolURL ' ' fullfile(toolDir, tool)]); end; addpath(genpath(fullfile(toolDir,tool)))

%%% matlab others
tool = 'red-blue-colormap'; repoURL = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/25536/versions/1/download/zip';
mathworksClone(repoURL, fullfile(toolDir, tool));
tool = 'util'; repoURL = 'https://github.com/Proulx-S/util.git'; subTool = ''; branch = 'dev-multiVencISMRM2026';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
tool = 'pcMRAsim'; repoURL = 'https://github.com/Proulx-S/pcMRAsim.git'; subTool = ''; branch = 'dev-multiVencISMRM2026';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
%% %%%%%%%%%%%%%%%%%%
disp(projectCode)
disp(projectStorage)
disp(projectScratch)
info.project.code            = projectCode;            clear projectCode
info.project.storage         = projectStorage;         clear projectStorage
info.project.scratch         = projectScratch;         clear projectScratch
info.project.dataBase        = projectDataBase;        clear projectDataBase
info.project.dataBasePhantom = projectDataBasePhantom; clear projectDataBasePhantom
info.toClean = {};





if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot simulation summary -- velocity map, mag map and and complex-domain signal evolution, for plug flow and laminar flow (both with flat magnitude profile) -- ISMRM2026-poster.pptx slide 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pSim7 = runSim; % get default parameters

% Flat magnitude profile: S.lumen = constant Mxy at vMean (no velocity-dependent inflow saturation)
Mz_flat7  = getMz_ss(pSim7.pMri, pSim7.pMri.relax.blood, pSim7.pVessel.vMean);
Mxy_flat7 = getMxy_ss(Mz_flat7, pSim7.pMri, pSim7.pMri.relax.blood);

pVesselPara7          = pSim7.pVessel;
pVesselPara7.S.lumen  = Mxy_flat7;

pVesselPlug7          = pSim7.pVessel;
pVesselPlug7.PD       = pVesselPlug7.ID;
pVesselPlug7.S.lumen  = Mxy_flat7;

% PCmono with venc list matching in vivo data range
pMri7 = pSim7.pMri;
pMri7.venc.method  = 'PCmono';
pMri7.venc.vencList = [40 20 13 10 8 7 6 5 4]';
pMri7.venc.FVEres = 0; pMri7.venc.FVEbw = 0;
pMri7.venc.FVEvel = []; pMri7.venc.vencMin = []; pMri7.venc.vencMax = [];

% Run simulations (light=false to retain magMap/vMap)
resPara7 = runSim(pVesselPara7, pSim7.pSim, pMri7, [], false);
resPlug7 = runSim(pVesselPlug7, pSim7.pSim, pMri7, [], false);

% Build multi-venc signal arrays for plotMultiVenc
% res.I dims: [1 1 1 1 nVenc 2] = [FE PE SL t M1 M1ref]
% col 1 = velocity-encoded, col 2 = M1=0 reference (real after phase subtraction)
vencListSim7 = pMri7.venc.vencList;
resListSim7  = {resPlug7, resPara7};
flowNamesSim7 = {'plug flow','laminar flow'};
IplotSim7    = cell(1,2);
vPlotSim7    = cell(1,2);
for flowIdx = 1:2
    r     = resListSim7{flowIdx};
    Iref  = r.I(1,1,1,1,1,2);
    Ienc  = squeeze(r.I(1,1,1,1,:,1));
    IplotSim7{flowIdx} = [Iref; Ienc(:)];
    vPlotSim7{flowIdx} = [inf;  vencListSim7];
end

% Figure: 2 rows (plug | laminar) x 3 cols (magMap | vMap | complex-plane spiral)
fSim7 = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 22]);
hTSim7 = tiledlayout(fSim7,2,3,'TileSpacing','compact','Padding','compact'); axSim7 = {};
for rowIdx = 1:2
    r_ = resListSim7{rowIdx};
    axSim7{end+1} = nexttile(hTSim7);
    imagesc(r_.magMap); axis image;
    axSim7{end}.Colormap = gray;
    axSim7{end}.CLim = [0 max(r_.magMap(:))*1.1];
    set(axSim7{end},'XTick',[],'YTick',[]);
    title(axSim7{end},{flowNamesSim7{rowIdx},'magnitude map'});
    axSim7{end+1} = nexttile(hTSim7);
    vLim7 = max(abs(r_.vMap(:)))*1.1;
    imagesc(r_.vMap,[-vLim7 vLim7]); axis image;
    axSim7{end}.Colormap = redblue;
    set(axSim7{end},'XTick',[],'YTick',[]);
    ylabel(colorbar,'velocity (cm/s)','FontSize',8);
    title(axSim7{end},{flowNamesSim7{rowIdx},'velocity map'});
    axSim7{end+1} = nexttile(hTSim7);
    runArr = ones(length(IplotSim7{rowIdx}),1);
    plotMultiVenc(axSim7{end}, IplotSim7{rowIdx}, vPlotSim7{rowIdx}, runArr, 'tight', 'jet');
    title(axSim7{end},{flowNamesSim7{rowIdx},'complex-domain signal'});
end
if ~exist(fullfile(info.project.storage, 'figures'),'dir'); mkdir(fullfile(info.project.storage, 'figures')); end
if saveThis || ~exist(fullfile(info.project.storage,'figures','simSummary.fig'),'file')
    saveas(        fSim7, fullfile(info.project.storage,'figures','simSummary.fig'));
    exportgraphics(fSim7, fullfile(info.project.storage,'figures','simSummary.png'));
    exportgraphics(fSim7, fullfile(info.project.storage,'figures','simSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




forceThis = 1;
%%%%%%%%%%%%%%%%%%%%
%% Load phantom data
%%%%%%%%%%%%%%%%%%%%
phantom03dataFile = fullfile(info.project.scratch, 'phantom03.mat');
if forceThis || ~exist(phantom03dataFile,'file')
    [data, dataVenc, dataRun, dataMeas, dataNoFlow, dataNoFlowMeas, PEspacing, FEspacing] = loadPhantom03(fullfile(info.project.dataBasePhantom,'20251010_multiVENCphantom03'));
    save(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
else
    load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
end

ID =  6.35; % mm
OD = 11.11; % mm

bestVenc = 10; % cm/s

% Flip phase sign
data = conj(data);
dataNoFlow = conj(dataNoFlow);

% % Compute PD and CD (meas by meas)
% runList = unique(dataRun);
% dataPD = nan(size(data));
% dataCD = nan(size(data));
% dataNoFlowPD = nan(size(dataNoFlow));
% dataNoFlowCD = nan(size(dataNoFlow));
% for rIdx = 1:length(runList)

%     measList = unique(dataMeas(dataRun==runList(rIdx)));
%     for mIdx = 1:length(measList)
%         idx    = dataRun==runList(rIdx) & dataMeas==measList(mIdx);
%         idxRef = dataRun==runList(rIdx) & dataMeas==measList(mIdx) & dataVenc==inf;
%         dataPD(:,:,idx)   = angle(  data(:,:,idx) ./ exp(1j.*angle(data(:,:,idxRef)))  );
%         dataCD(:,:,idx)   = data(:,:,idx) - data(:,:,idxRef);
%     end

%     measList = unique(dataNoFlowMeas(dataRun==runList(rIdx)));
%     for mIdx = 1:length(measList)
%         idx    = dataRun==runList(rIdx) & dataNoFlowMeas==measList(mIdx);
%         idxRef = dataRun==runList(rIdx) & dataNoFlowMeas==measList(mIdx) & dataVenc==inf;
%         dataNoFlowPD(:,:,idx)   = angle(  dataNoFlow(:,:,idx) ./ exp(1j.*angle(dataNoFlow(:,:,idxRef)))  );
%         dataNoFlowCD(:,:,idx)   = dataNoFlow(:,:,idx) - dataNoFlow(:,:,idxRef);
%     end
% end


% Compute coordinates around center of mass
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, size(M,1)*FEspacing-FEspacing/2, size(M,1));
PEpos = linspace(PEspacing/2, size(M,2)*PEspacing-PEspacing/2, size(M,2));
[FEgrid, PEgrid] = ndgrid(FEpos, PEpos);
total = sum(M(:));
com(1) = sum(FEgrid(:) .* M(:)) / total;  % center of mass along x (columns)
com(2) = sum(PEgrid(:) .* M(:)) / total;  % center of mass along y (rows)
FEgrid = FEgrid - com(1);
FEpos  = FEpos  - com(1);
PEgrid = PEgrid - com(2);
PEpos  = PEpos  - com(2);
rGrid = sqrt(PEgrid.^2+FEgrid.^2);
clear com


% Compute some masks
dFE = abs(FEgrid);
dPE = abs(PEgrid);
% farthest corner of each pixel from the center
d_far  = sqrt((dFE + FEspacing/2).^2 + (dPE + PEspacing/2).^2);
% nearest point of each pixel to the center
d_near = sqrt(max(0, dFE - FEspacing/2).^2 + max(0, dPE - PEspacing/2).^2);                                                                   

maskBloodOnly  = d_far  < ID/2;     % pixel entirely inside inner circle
maskTissueOnly = d_near > OD/2;    % pixel entirely outside outer circle  
maskWallLowMag = single(M<0.44e-7); % low magnitude pixels

theta = linspace(0, 2*pi, 360);                                                                                                               
%% %%%%%%%%%%%%%%%%%



if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom details -- all maps and masks with and without flow -- ISMRM2026-poster.pptx slide 8 and 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mag flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); M = {};
cLim = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    M{end+1} = abs(mean(data(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,M{end}); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        ylabel(colorbar(ax{end},'Location','westoutside'), 'MR magn. [a.u.]');
    end
    ax{end}.Colormap = gray;
    cLim{vencIdx} = clim(ax{end});
end
cLim = [0 max(max([M{:}]))];
set([ax{:}],'CLim',cLim);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'magFlowOn.fig'),'file')
    if ~exist(fullfile(info.project.storage, 'figures'),'dir'); mkdir(fullfile(info.project.storage, 'figures')); end
    saveas(        f      , fullfile(info.project.storage, 'figures', 'magFlowOn.fig'       ));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'magFlowOn.png'       ));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'magFlowOn.svg'       ));
end



% Mag flow off
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); M = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    M{end+1} = abs(mean(dataNoFlow(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,M{end}); axis image
    title(ax{end},['flow off; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        ylabel(colorbar(ax{end},'Location','westoutside'), 'MR magn. [a.u.]');
    end
    ax{end}.Colormap = gray;
end
set([ax{:}],'CLim',cLim);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');
    
%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'magFlowOff.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'magFlowOff.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'magFlowOff.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'magFlowOff.svg'));
end


% Phase flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); P = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    P{end+1} = angle(mean(data(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,P{end},[-pi pi]); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'PD [rad]');
    end
    ax{end}.Colormap = redblue;
    cb.Ticks = -pi:pi/2:pi;
    cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'}; 
end

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'phaseFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'phaseFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phaseFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phaseFlowOn.svg'));
end



% Phase flow off
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); P = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    P{end+1} = angle(mean(dataNoFlow(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,P{end},[-pi pi]); axis image
    title(ax{end},['flow off; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'PD [rad]');
    end
    ax{end}.Colormap = redblue;
    cb.Ticks = -pi:pi/2:pi;
    cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'}; 
end

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'phaseFlowOff.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'phaseFlowOff.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phaseFlowOff.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phaseFlowOff.svg'));
end




% CDvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); PHI = {}; CDvel = {};
cLim = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    % P{end+1} = angle(mean(data(:,:,dataVenc==vencList(vencIdx)),3));
    CD = mean(data(:,:,dataVenc==vencList(vencIdx)),3) - mean(data(:,:,dataVenc==inf),3);
    [CDvel{end+1},PHI{end+1}] = getPlugFlowEstimates(vencList(vencIdx),CD,[],[],[],0);
    CDvel{end}(maskWallLowMag | rGrid>(ID/2+OD/2)./2) = nan;
    
    
    % imagesc(ax{end},PEpos,FEpos,PHI{end},[-pi pi]); axis image
    imagesc(ax{end},PEpos,FEpos,CDvel{end}); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        % ylabel(cb, 'phi [rad]');
        ylabel(cb, 'velocity [cm/s]');
    end
    ax{end}.Colormap = redblue;
    % cb.Ticks = -pi:pi/2:pi;
    % cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'}; 
    cLim{vencIdx} = clim(ax{end});
end
% cLim = abs([CDvel{:}]); cLim = [-1 1].*max((tmp(cLim~=inf)));
set([ax{:}],'CLim',[-9 9]);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'CDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'CDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'CDvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'CDvelFlowOn.svg'));
end



% PDvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); PDvel = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    % P{end+1} = angle(mean(data(:,:,dataVenc==vencList(vencIdx)),3));
    PD = angle(mean(  data(:,:,dataVenc==vencList(vencIdx)) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3)))  ,3));
    PDvel{end+1} = phase2vel(PD,vencToM1(vencList(vencIdx)));
    imagesc(ax{end},PEpos,FEpos,PDvel{end}); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'PD velocity [cm/s]');
    end
    ax{end}.Colormap = redblue;
    % cb.Ticks = -pi:pi/2:pi;
    % cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'}; 
end
set([ax{:}],'CLim',[-1 1].*9);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'PDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'PDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'PDvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'PDvelFlowOn.svg'));
end




% Pvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); Pvel = {};
for vencIdx = 1:size(vencList,1)
    ax{end+1} = nexttile;
    P = angle(mean(  data(:,:,dataVenc==vencList(vencIdx))  ,3));
    Pvel{end+1} = phase2vel(P,vencToM1(vencList(vencIdx)));
    imagesc(ax{end},PEpos,FEpos,Pvel{end}); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'PD velocity [cm/s]');
    end
    ax{end}.Colormap = redblue;
    % cb.Ticks = -pi:pi/2:pi;
    % cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'}; 
end
set([ax{:}],'CLim',[-1 1].*9);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');

%save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'PvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.storage, 'figures', 'PvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'PvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'PvelFlowOn.svg'));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom details -- radial profiles -- ISMRM2026-poster.pptx slide 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_rad  = abs(mean(data(:,:,dataVenc==inf),3));
PD_rad = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))), 3));
m1best = vencToM1(bestVenc);
velPD_rad = phase2vel(PD_rad, m1best); % [cm/s]
% Dominant flow sign after conj() may be negative; flip so positive = into slice
if mean(velPD_rad(maskBloodOnly)) < 0; velPD_rad = -velPD_rad; end

idxBlood = maskBloodOnly;
idxSel1  = maskBloodOnly | (maskTissueOnly & rGrid<OD/2);

% Parabolic fit to velocity vs. radius (blood-only pixels)
rFit3  = double(rGrid(idxBlood));
vFit3  = double(velPD_rad(idxBlood));
parabola_fun3 = @(vMax, r) vMax .* (1 - (r./(ID/2)).^2);
rFine3 = linspace(0, ID/2, 100);
lb3 = 0; opts3 = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000);
vMax_fit3 = lsqcurvefit(parabola_fun3, max(vFit3), rFit3(:), vFit3(:), lb3, [], opts3);
vMean_fit3 = vMax_fit3 / 2;
vFitLine3  = parabola_fun3(vMax_fit3, rFine3);
fprintf('Radial profile fit: vMax=%.3f cm/s, vMean=%.3f cm/s\n', vMax_fit3, vMean_fit3);

% Polynomial fit to magnitude vs. radius (blood-only, linear in r^2)
pMag3     = polyfit(double(rGrid(idxBlood)).^2, double(M_rad(idxBlood)), 1);
MFitLine3 = polyval(pMag3, rFine3.^2);

theta_circ = linspace(0,2*pi,360);

f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact'); ax = {};

% velocity vs. radial position with parabolic fit
ax{end+1} = nexttile(hT);
plot(rGrid(idxSel1), velPD_rad(idxSel1), '.', 'Color', [0.5 0.5 0.8]);
hold on
plot(rFine3, vFitLine3, 'r-', 'LineWidth', 2);
xlabel('off-center position [mm]'); ylabel('velocity (cm/s)');
title(sprintf('velocity profile\nvMax=%.2f, vMean=%.2f cm/s', vMax_fit3, vMean_fit3));
grid on; legend('data','parabolic fit','Location','north');

% magnitude vs. radial position with polynomial fit
ax{end+1} = nexttile(hT);
plot(rGrid(idxSel1), M_rad(idxSel1), '.', 'Color', [0.5 0.8 0.5]);
hold on
plot(rFine3, MFitLine3, 'r-', 'LineWidth', 2);
xlabel('off-center position [mm]'); ylabel('MR signal magnitude [a.u.]');
title('magnitude profile'); grid on;
legend('data','polynomial fit','Location','north');

% phase vs. magnitude scatter
ax{end+1} = nexttile(hT);
plot(velPD_rad(idxSel1), M_rad(idxSel1), '.', 'Color', [0.8 0.5 0.5]);
xlabel('velocity (cm/s)'); ylabel('MR signal magnitude [a.u.]');
title('velocity vs. magnitude'); grid on;

% velocity 2D map
ax{end+1} = nexttile(hT);
vLim3 = max(abs(velPD_rad(:)));
imagesc(PEpos, FEpos, velPD_rad, [-vLim3 vLim3]); axis image;
ax{end}.Colormap = redblue;
hold on; plot(ID/2*cos(theta_circ), ID/2*sin(theta_circ), 'w--', 'LineWidth', 1);
set(ax{end},'XTick',[],'YTick',[]); title(sprintf('velocity map\n(venc=%gcm/s)',bestVenc));
ylabel(colorbar,'velocity [cm/s]');

% magnitude 2D map
ax{end+1} = nexttile(hT);
imagesc(PEpos, FEpos, M_rad, [0 max(M_rad(:))]); axis image;
ax{end}.Colormap = gray;
hold on; plot(ID/2*cos(theta_circ), ID/2*sin(theta_circ), 'w--', 'LineWidth', 1);
set(ax{end},'XTick',[],'YTick',[]); title('magnitude map (venc=inf)');

% blood-only mask
ax{end+1} = nexttile(hT);
imagesc(PEpos, FEpos, maskBloodOnly, [0 1]); axis image;
ax{end}.Colormap = gray;
hold on; plot(ID/2*cos(theta_circ), ID/2*sin(theta_circ), 'm', 'LineWidth', 1.5);
set(ax{end},'XTick',[],'YTick',[]); title('blood-only mask');

if saveThis || ~exist(fullfile(info.project.storage,'figures','radialProfiles.fig'),'file')
    saveas(        f, fullfile(info.project.storage,'figures','radialProfiles.fig'));
    exportgraphics(f, fullfile(info.project.storage,'figures','radialProfiles.png'));
    exportgraphics(f, fullfile(info.project.storage,'figures','radialProfiles.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom summary -- reference mag, velocity map at good venc and complex-domain signal evolution -- ISMRM2026-poster.pptx slide 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT = tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};

% Plot mag
ax{end+1} = nexttile;
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
hIm = imagesc(ax{end},PEpos,FEpos,M,[0 max(M(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray;
set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},'phantom ROI');

% Plot velocity map
ax{end+1} = nexttile;
PD   = angle(mean(  data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3)))  ,3));
CD   = mean(data(:,:,dataVenc==bestVenc),3)-mean(data(:,:,dataVenc==inf),3);
[velCD,phi,velPD,~,~,~] = getPlugFlowEstimates(bestVenc,CD,[],[],PD);
hIm = imagesc(ax{end},PEpos,FEpos,velPD,[-max(abs(velPD(:))) max(abs(velPD(:)))]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = redblue;
set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},['venc=' num2str(bestVenc) 'cm/s']);

% Plot masks
ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly,[0 1]); axis image
ax{end}.Colormap = gray;
set(ax{end},'XTick',[],'YTick',[],'Color','none');
hold(ax{end},'on');
theta = linspace(0, 2*pi, 360);
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
title(ax{end},'blood-only mask');
ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly,[0 1]); axis image
ax{end}.Colormap = gray;
set(ax{end},'XTick',[],'YTick',[],'Color','none');
hold(ax{end},'on');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'tissue-only mask');
ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskWallLowMag,[0 1]); axis image
ax{end}.Colormap = gray;
set(ax{end},'XTick',[],'YTick',[],'Color','none');
hold(ax{end},'on');
plot(ax{end},ID/2 * cos(theta), ID/2 * sin(theta), 'm');
plot(ax{end},OD/2 * cos(theta), OD/2 * sin(theta), 'm');
title(ax{end},'low-mag wall mask');

% Plot complex domain signal
ax{end+1} = nexttile([3 3]);
I     = squeeze(mean(data,[1 2]));
Ivenc = squeeze(dataVenc);
Irun  = squeeze(dataRun);
plotMultiVenc(ax{end},I,Ivenc,Irun,[],'hot');

% Save
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'phantomSummary.fig'))
    saveas(        f      , fullfile(info.project.storage, 'figures', 'phantomSummary.fig'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phantomSummary.png'));
    exportgraphics(f      , fullfile(info.project.storage, 'figures', 'phantomSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot matched simulation summary -- same format as phantom summary but all from simulation
%%   -- spatial maps (rows=PE, cols=FE in runSim → transposed for display: rows=FE, cols=PE)
%%   -- magnitude calibrated so lumen/surround ratio matches phantom
%%   -- complex-domain spiral with fine, evenly-M1-spaced venc sweep (smooth line)
%%   -- ISMRM2026-poster.pptx slide 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine M1-spaced venc list for smooth spiral (400 steps, evenly spaced in M1)
gamma_phys_sim = 2.6752218708e8 / (2*pi); % Hz/T
M1_max_sim     = vencToM1(2);             % M1 for venc=2 cm/s
nSteps_sim     = 400;
M1_fine_sim    = linspace(M1_max_sim/nSteps_sim, M1_max_sim, nSteps_sim);
venc_fine_sim  = pi * 100 ./ (gamma_phys_sim .* M1_fine_sim); % [cm/s]

pDefSim = runSim;
pVesselSim           = pDefSim.pVessel;
pVesselSim.ID        = ID;
pVesselSim.WT        = 2.38125;
pVesselSim.vMean     = vMean_fit3;
pVesselSim.profile   = 'parabolic1';
pVesselSim.S.lumen   = []; % velocity-dependent inflow (auto)
pVesselSim.S.wall    = 0;
pSimSim              = pDefSim.pSim;
pSimSim.fovFE        = size(data,1) * FEspacing; % match phantom FOV exactly
pSimSim.fovPE        = size(data,2) * PEspacing;
pSimSim.matFE        = 3;
pSimSim.matPE        = 3;
pSimSim.nSpin        = (2^8+1)^2;
pMriSim              = pDefSim.pMri;
pMriSim.venc.method  = 'PCmono';
pMriSim.venc.vencList = venc_fine_sim(:);
pMriSim.venc.FVEres  = 0; pMriSim.venc.FVEbw = 0;
pMriSim.venc.FVEvel  = []; pMriSim.venc.vencMin = []; pMriSim.venc.vencMax = [];
resSim = runSim(pVesselSim, pSimSim, pMriSim, [], false); % light=false → keep magMap/vMap

% Calibrate surround magnitude to match phantom lumen/tissue ratio
lumen_mask_sim    = resSim.pVessel.mask.lumen;
wall_mask_sim     = resSim.pVessel.mask.wall;
surround_mask_sim = resSim.pVessel.mask.surround;
magMap_sim = double(resSim.magMap); % per-spin units (S/nSpin)
S_lumen_ps = mean(magMap_sim(lumen_mask_sim));
S_surround_corrected_ps = S_lumen_ps / (mean(abs(M(maskBloodOnly))) / mean(abs(M(maskTissueOnly))));
magMap_sim(surround_mask_sim) = S_surround_corrected_ps;
magMap_sim(wall_mask_sim)     = 0;

% Grid axes: runSim grid rows=PE, cols=FE → transpose for display (rows=FE, cols=PE)
FEax_sim = resSim.pSim.gridFE(1,:); % FE values (varies along cols)
PEax_sim = resSim.pSim.gridPE(:,1); % PE values (varies along rows)
magMap_T    = magMap_sim';           % after transpose: rows=FE, cols=PE
vMap_T      = double(resSim.vMap)';
lumen_T     = lumen_mask_sim';
wall_T      = wall_mask_sim';
surround_T  = surround_mask_sim';

% Complex-plane spiral signal
Iref_sim      = resSim.I(1,1,1,1,1,2);
Ienc_sim      = squeeze(resSim.I(1,1,1,1,:,1));
Ienc_norm_sim = Ienc_sim / abs(Iref_sim);
theta_sim     = linspace(0,2*pi,360);

fSim = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hTSim = tiledlayout(fSim,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); axSim = {};

% Magnitude map
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, magMap_T, [0 max(magMap_T(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end},'simulation ROI');

% Velocity map
axSim{end+1} = nexttile(hTSim);
vLim_sim = max(abs(vMap_T(:)));
imagesc(axSim{end}, PEax_sim, FEax_sim, vMap_T, [-vLim_sim vLim_sim]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axSim{end}.Colormap = redblue; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end}, sprintf('parabolic vMean=%.1f cm/s', vMean_fit3));

% Lumen mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(lumen_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on'); plot(axSim{end}, ID/2*cos(theta_sim), ID/2*sin(theta_sim), 'm');
title(axSim{end},'lumen mask');

% Wall mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(wall_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on');
plot(axSim{end}, ID/2*cos(theta_sim),       ID/2*sin(theta_sim), 'm');
plot(axSim{end}, (ID/2+2.38125)*cos(theta_sim), (ID/2+2.38125)*sin(theta_sim), 'm');
title(axSim{end},'wall mask');

% Surround mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(surround_T), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on'); plot(axSim{end}, (ID/2+2.38125)*cos(theta_sim), (ID/2+2.38125)*sin(theta_sim), 'm');
title(axSim{end},'surround mask');

% Complex-plane spiral: smooth colored line (surface trick for per-segment coloring)
axSim{end+1} = nexttile(hTSim, [3 3]);
xSp  = real(Ienc_norm_sim);
ySp  = imag(Ienc_norm_sim);
cVal = (1:nSteps_sim)' / nSteps_sim;
surface(axSim{end}, [xSp xSp], [ySp ySp], zeros(nSteps_sim,2), [cVal cVal], ...
    'EdgeColor','flat','FaceColor','none','LineWidth',1.5);
colormap(axSim{end}, jet);
cb_sim = colorbar(axSim{end},'Location','eastoutside');
ylabel(cb_sim,'normalized M_1  (0=low → 1=high)');
hold(axSim{end},'on');
plot(axSim{end}, cos(theta_sim), sin(theta_sim), 'w--', 'LineWidth', 0.8);
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

if saveThis || ~exist(fullfile(info.project.storage,'figures','matchedSimSummary.fig'),'file')
    saveas(        fSim, fullfile(info.project.storage,'figures','matchedSimSummary.fig'));
    exportgraphics(fSim, fullfile(info.project.storage,'figures','matchedSimSummary.png'));
    exportgraphics(fSim, fullfile(info.project.storage,'figures','matchedSimSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in vivo data -- sub-01 and sub-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI coordinates (vessels defined in multiVencInVivo project)
inVivoSubRoiList = {};
% sub-01: 6 vessels
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY = [37 47];   inVivoSubRoiList{end}(1).roiX = [87 92];
inVivoSubRoiList{end}(2).roiY = [158 164]; inVivoSubRoiList{end}(2).roiX = [90 92];
inVivoSubRoiList{end}(3).roiY = [91 94];   inVivoSubRoiList{end}(3).roiX = [88 91];
inVivoSubRoiList{end}(4).roiY = [103 107]; inVivoSubRoiList{end}(4).roiX = [77 81];
inVivoSubRoiList{end}(5).roiY = [100 103]; inVivoSubRoiList{end}(5).roiX = [139 141];
inVivoSubRoiList{end}(6).roiY = [130 134]; inVivoSubRoiList{end}(6).roiX = [50 53];
% sub-02: 5 vessels
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY = [235 242]; inVivoSubRoiList{end}(1).roiX = [140 144];
inVivoSubRoiList{end}(2).roiY = [228 232]; inVivoSubRoiList{end}(2).roiX = [55 59];
inVivoSubRoiList{end}(3).roiY = [199 203]; inVivoSubRoiList{end}(3).roiX = [136 140];
inVivoSubRoiList{end}(4).roiY = [216 219]; inVivoSubRoiList{end}(4).roiX = [52 54];
inVivoSubRoiList{end}(5).roiY = [163 169]; inVivoSubRoiList{end}(5).roiX = [91 94];

inVivoSubNames   = {'sub-01','sub-02'};
inVivoScratch    = fullfile(fileparts(info.project.code), 'multiVencInVivo', 'tmp');
inVivoSubData    = cell(1,2);
for s = 1:2
    subFile = fullfile(inVivoScratch, [inVivoSubNames{s} '.mat']);
    inVivoSubData{s} = load(subFile, 'img', 'imgInfo', 'refImgAv');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot in vivo summary -- reference mag, velocity map and complex-domain signal evolution -- ISMRM2026-poster.pptx slide 10 (sub-01) and 11 (sub-02)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fullfile(info.project.storage, 'figures'),'dir'); mkdir(fullfile(info.project.storage, 'figures')); end

for s = 1:2
    img      = inVivoSubData{s}.img;
    imgInfo  = inVivoSubData{s}.imgInfo;
    refImgAv = inVivoSubData{s}.refImgAv;

    % ROI overlay figure
    hF = figure('MenuBar','none','ToolBar','none','Units','normalized','Position',[0 0 1 1]);
    hIm = imagesc(mean(abs(refImgAv),7));
    for roiIdx = 1:length(inVivoSubRoiList{s})
        roiY = inVivoSubRoiList{s}(roiIdx).roiY;
        roiX = inVivoSubRoiList{s}(roiIdx).roiX;
        hold on
        hBox(roiIdx)  = plot([roiX(1)-0.5 roiX(2)+0.5 roiX(2)+0.5 roiX(1)-0.5 roiX(1)-0.5], ...
                             [roiY(1)-0.5 roiY(1)-0.5 roiY(2)+0.5 roiY(2)+0.5 roiY(1)-0.5], 'c');
        hText(roiIdx) = text(roiX(1),roiY(1),sprintf('roi%d',roiIdx),'Color','r','FontSize',12,'FontWeight','bold');
    end
    ax = gca; axis image;
    set(ax,'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
    ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    ax.CLim = [0, 1/2*max(hIm.CData(:))];
    set(hBox,'LineWidth',0.5);
    if saveThis || ~exist(fullfile(info.project.storage,'figures',[inVivoSubNames{s} '-roiOverlay.png']),'file')
        drawnow;
        exportgraphics(hF, fullfile(info.project.storage,'figures',[inVivoSubNames{s} '-roiOverlay.png']));
        exportgraphics(hF, fullfile(info.project.storage,'figures',[inVivoSubNames{s} '-roiOverlay.svg']));
    end
    close(hF); clear hBox hText

    % Per-vessel spiral figures
    for roiIdx = 1:length(inVivoSubRoiList{s})
        roiY = inVivoSubRoiList{s}(roiIdx).roiY;
        roiX = inVivoSubRoiList{s}(roiIdx).roiX;

        trj = img(roiY(1):roiY(2), roiX(1):roiX(2), :,:,:,:,:,:,:,:,:,:,:,:,:,:);

        runIdxList = unique(imgInfo.runIdx);
        for runIdx = 1:length(runIdxList)
            idx = squeeze(imgInfo.runIdx==runIdxList(runIdx) & imgInfo.vencList==inf);
            refPhase = angle(mean(trj(:,:,:,:,:,:,idx,:,:,:,:,:,:,:,:,:), [7 11]));
            idx2 = squeeze(imgInfo.runIdx==runIdxList(runIdx));
            trj(:,:,:,:,:,:,idx2,:,:,:,:,:,:,:,:,:) = ...
                trj(:,:,:,:,:,:,idx2,:,:,:,:,:,:,:,:,:) ./ exp(1i*refPhase);
        end

        trj     = permute(mean(trj,[1 2]),[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
        trj     = trj ./ abs(mean(trj(1:2,:),[1 2]));
        trjVenc = permute(imgInfo.vencList,[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
        trjVencLabel = replace(cellstr(num2str(trjVenc)),' ','');

        d = {};
        for vencIdx = 1:size(trj,1)
            d{vencIdx,1} = [real(trj(vencIdx,:)); imag(trj(vencIdx,:))]';
        end
        prob = credibleMean2d(d);

        hFv = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 20 20]);
        cMap_venc = sort(unique(trjVenc),'descend');
        cMap = jet(length(cMap_venc));
        for dIdx = 1:length(d)
            cMap_idx = cMap_venc==trjVenc(dIdx);
            hPcont(dIdx) = plot(prob{dIdx}.CIcontour,'FaceColor',cMap(cMap_idx,:),'EdgeColor','k');
            hold on
        end
        axv = gca;
        set(hPcont,'FaceAlpha',1);
        axis image tight
        xLim = xlim; if xLim(1)>0; xLim(1)=0; end; if xLim(2)<0; xLim(2)=0; end;
        yLim = ylim; if yLim(1)>0; yLim(1)=0; end; if yLim(2)<0; yLim(2)=0; end;
        dLim = max(diff(xLim),diff(yLim)).*0.03;
        xLim = xLim+[-dLim dLim]; yLim = yLim+[-dLim dLim];
        set(axv,'XLim',xLim,'YLim',yLim);
        uistack(xline(0,'w'),'bottom'); uistack(yline(0,'w'),'bottom');
        grid on; xlabel('Re','FontName','Times New Roman'); ylabel('Im','FontName','Times New Roman');
        axv.Color = 'k';
        dTick = min(mean(diff(axv.XTick)),mean(diff(axv.YTick)));
        xTick = round(xLim./dTick).*dTick; xTick = linspace(xTick(1),xTick(2),range(xTick)/dTick+1);
        yTick = round(yLim./dTick).*dTick; yTick = linspace(yTick(1),yTick(2),range(yTick)/dTick+1);
        set(axv,'XTick',xTick,'YTick',yTick);
        theta = linspace(0,2*pi,100);
        x = abs(mean(trj(trjVenc==inf,:),[1 2]))*cos(theta);
        y = abs(mean(trj(trjVenc==inf,:),[1 2]))*sin(theta);
        uistack(plot(x,y,'w'),'bottom');
        title(legend(hPcont,trjVencLabel),'V_{enc} (cm/s)');
        vesselFigName = [inVivoSubNames{s} '_vessel-' num2str(roiIdx,'%02d')];
        if saveThis || ~exist(fullfile(info.project.storage,'figures',[vesselFigName '.png']),'file')
            drawnow;
            exportgraphics(hFv, fullfile(info.project.storage,'figures',[vesselFigName '.png']));
            exportgraphics(hFv, fullfile(info.project.storage,'figures',[vesselFigName '.svg']));
        end
        close(hFv); clear hPcont d prob
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




if 1
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Illustrate the effect of inflow on FVE spectrum -- ISMRM2026-poster.pptx slide 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = runSim;
tmp = runSim(p.pVessel,p.pSim,p.pMri);
p.pVessel.vMean = tmp.pMri.vCrit/2*1.5;
p.pMri.FA = 90;
p.pMri.venc.method = 'FVEmono';
p.pMri.venc.FVEbw = p.pVessel.vMean*4;
p.pMri.venc.FVEres = p.pMri.venc.FVEbw./200;
p.pSim.nSpin = (2^10+1)^2;
res      = runSim(p.pVessel,p.pSim,p.pMri,[],0);
% p.pMri.sliceThickness = inf;
Mz  = getMz_ss(p.pMri,p.pMri.relax.blood,p.pVessel.vMean);
Mxy = getMxy_ss(Mz,p.pMri,p.pMri.relax.blood);
p.pVessel.S.lumen = Mxy;
resSatin = runSim(p.pVessel,p.pSim,p.pMri);

fVelSpec = figure;
hFlat = plot(resSatin.pMri.venc.FVEvel,abs(fftshift(fft(squeeze(resSatin.I)))),'w');
hold on
hVdep = plot(res.pMri.venc.FVEvel     ,abs(fftshift(fft(squeeze(res.I))))     ,'g');
axis tight; grid on; xlabel('velocity (cm/s)'); ylabel('velocity spectrum/histogram');
yLim = ylim; yLim(1) = 0; ylim(yLim);
[N,edges] = histcounts(res.vMap(res.pSim.gridVoxIdx==0),20);
% binWidth = mean(diff(edges));
% edges = edges-binWidth/2; edges(end+1) = edges(end)+binWidth;
% [N,edges] = histcounts(res.vMap(res.pSim.gridVoxIdx==0),edges);
hVhist = histogram('BinEdges',edges,'BinCounts',N/max(N)*yLim(2),'FaceColor',0.5.*[1 1 1],'EdgeColor','none');
legend([hFlat,hVdep,hVhist],['velocity spectrum from' newline 'flat magnitude profile'],['velocity spectrum from' newline 'velocity-dependent magnitude profile'],'normalized velocity histogram','Location','northwest','box','off');
uistack(hVhist,'bottom');

fVelSpecInflow = figure;
% inflowVel = linspace(0,p.pVessel.vMean*3,2^10);
inflowVel = linspace(0,6,2^10);
[inflowMz,~,~,~,~,inflowVel] = getMz_ss(p.pMri,p.pMri.relax.blood,inflowVel);
% inflowMxy = getMxy_ss(inflowMz,p.pMri,p.pMri.relax.blood);
hStairs = stairs(inflowVel,inflowMz,'g');
axis tight square; grid on; xlabel('spin velocity (cm/s)'); ylabel('M_z');
ylim([0 1])

if ~exist(fullfile(info.project.storage, 'figures'),'dir'); mkdir(fullfile(info.project.storage, 'figures')); end
if saveThis || ~exist(fullfile(info.project.storage, 'figures', 'FVEvelSpec.fig'),'file') || ~exist(fullfile(info.project.storage, 'figures', 'FVEvelSpec_inflow.fig'),'file')
    saveas(        fVelSpec      , fullfile(info.project.storage, 'figures', 'FVEvelSpec.fig'       ));
    exportgraphics(fVelSpec      , fullfile(info.project.storage, 'figures', 'FVEvelSpec.png'       ));
    exportgraphics(fVelSpec      , fullfile(info.project.storage, 'figures', 'FVEvelSpec.svg'       ));
    saveas(        fVelSpecInflow, fullfile(info.project.storage, 'figures', 'FVEvelSpec_inflow.fig'));
    exportgraphics(fVelSpecInflow, fullfile(info.project.storage, 'figures', 'FVEvelSpec_inflow.png'));
    exportgraphics(fVelSpecInflow, fullfile(info.project.storage, 'figures', 'FVEvelSpec_inflow.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVelSpec; % FVE spectra reflects spin velocity distribution, but weighted by velocity-dependent spin magnitude
fVelSpecInflow; % Here the weighting effect was maximized using a 90 flip angle for a linear magnitude function of velocity
end







