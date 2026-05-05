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
info.project.figures = fullfile(info.project.code, 'figures'); if ~exist(info.project.figures,'dir'); mkdir(info.project.figures); end
info.toClean = {};





if 0
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
    I_s7    = IplotSim7{rowIdx};
    v_s7    = vPlotSim7{rowIdx};
    Mnorm_s7 = abs(mean(I_s7(v_s7==inf)));
    plotComplexDomain(axSim7{end}, I_s7/Mnorm_s7, 'tight', 'line');
    title(axSim7{end},{flowNamesSim7{rowIdx},'complex-domain signal'});
end
if ~exist(info.project.figures,'dir'); mkdir(info.project.figures); end
if saveThis || ~exist(fullfile(info.project.figures,'simSummary.fig'),'file')
    saveas(        fSim7, fullfile(info.project.figures,'simSummary.fig'));
    exportgraphics(fSim7, fullfile(info.project.figures,'simSummary.png'));
    exportgraphics(fSim7, fullfile(info.project.figures,'simSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




forceThis = 0;
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

maskBloodOnly  = d_far  < ID/2;              % pixel entirely inside inner circle
maskWallOnly   = d_near > ID/2 & d_far < OD/2; % pixel entirely within wall annulus
maskTissueOnly = d_near > OD/2;              % pixel entirely outside outer circle
maskWallLowMag = single(M<0.44e-7);          % low magnitude pixels

theta = linspace(0, 2*pi, 360);                                                                                                               
%% %%%%%%%%%%%%%%%%%



if 0
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
if saveThis || ~exist(fullfile(info.project.figures,'magFlowOn.fig'),'file')
    if ~exist(info.project.figures,'dir'); mkdir(info.project.figures); end
    saveas(        f      , fullfile(info.project.figures,'magFlowOn.fig'       ));
    exportgraphics(f      , fullfile(info.project.figures,'magFlowOn.png'       ));
    exportgraphics(f      , fullfile(info.project.figures,'magFlowOn.svg'       ));
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
if saveThis || ~exist(fullfile(info.project.figures,'magFlowOff.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'magFlowOff.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'magFlowOff.png'));
    exportgraphics(f      , fullfile(info.project.figures,'magFlowOff.svg'));
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
if saveThis || ~exist(fullfile(info.project.figures,'phaseFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'phaseFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'phaseFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.figures,'phaseFlowOn.svg'));
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
if saveThis || ~exist(fullfile(info.project.figures,'phaseFlowOff.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'phaseFlowOff.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'phaseFlowOff.png'));
    exportgraphics(f      , fullfile(info.project.figures,'phaseFlowOff.svg'));
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
if saveThis || ~exist(fullfile(info.project.figures,'CDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'CDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'CDvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.figures,'CDvelFlowOn.svg'));
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
if saveThis || ~exist(fullfile(info.project.figures,'PDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'PDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'PDvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.figures,'PDvelFlowOn.svg'));
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
if saveThis || ~exist(fullfile(info.project.figures,'PvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(info.project.figures,'PvelFlowOn.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'PvelFlowOn.png'));
    exportgraphics(f      , fullfile(info.project.figures,'PvelFlowOn.svg'));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom details -- radial profiles -- ISMRM2026-poster.pptx slide 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_rad  = abs(mean(data(:,:,dataVenc==inf),3));
M_rad2 = abs(mean(data(:,:,dataVenc==2  ),3));
PD_rad = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))), 3));
m1best = vencToM1(bestVenc);
velPD_rad = phase2vel(PD_rad, m1best); % [cm/s]
% Dominant flow sign after conj() may be negative; flip so positive = into slice
if mean(velPD_rad(maskBloodOnly)) < 0; velPD_rad = -velPD_rad; end

idxBlood = maskBloodOnly;
idxVel   = ~logical(maskWallLowMag);
idxMag   = true(size(rGrid));

% Parabolic fit to velocity vs. radius (blood-only pixels)
rFit3  = double(rGrid(idxBlood));
vFit3  = double(velPD_rad(idxBlood));
rFine3 = linspace(0, max(rFit3), 100);
ft_vel  = fittype('vMax * (1 - (x/R)^2)', 'independent', 'x', 'problem', 'R', 'coefficients', {'vMax'});
vel_fit = fit(rFit3(:), vFit3(:), ft_vel, 'problem', {ID/2}, 'Lower', 0, 'StartPoint', max(vFit3));
vMax_fit3  = vel_fit.vMax;
vMean_fit3 = vMax_fit3 / 2;
vFitLine3  = vel_fit(rFine3(:));

ft_vel2  = fittype('vMax * (1 - (x/R)^2)', 'independent', 'x', 'coefficients', {'vMax', 'R'});
vel_fit2 = fit(rFit3(:), vFit3(:), ft_vel2, 'Lower', [0, 0], 'StartPoint', [max(vFit3), ID/2]);
vMax_fit3b  = vel_fit2.vMax;
R_fit3b     = vel_fit2.R;
vMean_fit3b = vMax_fit3b / 2;
vFitLine3b  = vel_fit2(rFine3(:));

fprintf('Constrained fit (R=ID/2=%.3fmm): vMax=%.3f cm/s, vMean=%.3f cm/s\n', ID/2, vMax_fit3, vMean_fit3);
fprintf('Free fit        (R=%.3fmm):      vMax=%.3f cm/s, vMean=%.3f cm/s\n', R_fit3b, vMax_fit3b, vMean_fit3b);

% % Polynomial fit to magnitude vs. radius (blood-only, linear in r^2)
% mag_fit   = fit(double(rGrid(idxBlood(:))).^2, double(M_rad(idxBlood(:))), 'poly1');
% MFitLine3 = mag_fit(rFine3(:).^2);

% Degree-2 polynomial in r with maximum fixed at r=0: f(r) = a + b*r^2
ft_mag  = fittype('a + b*x^2', 'independent', 'x', 'coefficients', {'a', 'b'});
mag_fit = fit(double(rGrid(idxBlood(:))), double(M_rad(idxBlood(:))), ft_mag, ...
    'StartPoint', [max(double(M_rad(idxBlood(:)))), -1]);
MFitLine3 = mag_fit(rFine3(:));

% Theoretical inflow enhancement: Mz as a function of velocity
p_inflow = runSim;
p_inflow.pMri.fieldStrength = 3;
p_inflow.pMri.species = 'phantom';
p_inflow.pMri.sliceThickness = 2.2;
p_inflow.pMri.TR = 75.90/(5+1)/1000;
p_inflow.pMri.TE = 9.8/1000;
p_inflow.pMri.FA = 50;
p_inflow = runSim([],[],p_inflow.pMri);

[Mz_inflow,~,~,~,~,velInflow] = getMz_ss(p_inflow.pMri, p_inflow.pMri.relax.blood);

theta_circ = linspace(0,2*pi,360);

f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact'); ax = {};

% velocity vs. radial position with parabolic fit
ax{end+1} = nexttile(hT);
plot(rGrid(idxVel), velPD_rad(idxVel), '.', 'Color', [0.5 0.5 0.8]);
hold on
% plot(rFine3, vFitLine3,  '-', 'Color', [0.5 0.5 0.8], 'LineWidth', 2);
plot(rFine3, vFitLine3b, '-', 'Color', [0.9 0.6 0.1], 'LineWidth', 2);
xline(ID/2, 'w--'); xline(OD/2, 'w--');
xlabel('off-center position [mm]'); ylabel('velocity (cm/s)');
title(sprintf('velocity profile\nvMean: est=%.2f, true=%.2f cm/s\ndiam: est=%.2f, true=%.2f mm', vMean_fit3b, vMean_fit3, 2*R_fit3b, ID));
grid on; legend('data', sprintf('free fit (R=%.2fmm)', R_fit3b), 'ID/2','OD/2','Location','north');

% magnitude vs. radial position with polynomial fit
ax{end+1} = nexttile(hT);
plot(rGrid(idxMag), M_rad( idxMag), '.', 'Color', [0.5 0.8 0.5]);
hold on
plot(rGrid(idxMag), M_rad2(idxMag), '.', 'Color', [0.8 0.6 0.2]);
plot(rFine3, MFitLine3, '-', 'Color', [0.5 0.8 0.5], 'LineWidth', 2);
xline(ID/2, 'w--'); xline(OD/2, 'w--');
xlabel('off-center position [mm]'); ylabel('MR signal magnitude [a.u.]');
title('magnitude profile'); grid on;
legend('venc=inf','venc=2','polynomial fit (venc=inf)','ID/2','OD/2','Location','north');

% velocity vs. magnitude scatter + theoretical inflow enhancement
ax{end+1} = nexttile(hT);
% yyaxis left
plot(velPD_rad(idxBlood), M_rad(idxBlood), '.', 'Color', [0.8 0.5 0.5]); hold on
xlabel('velocity (cm/s)'); ylabel('MR signal magnitude [a.u.]');
% ylim([0 inf]);
% yyaxis right
Mxy_inflow = getMxy_ss(Mz_inflow, p_inflow.pMri, p_inflow.pMri.relax);
stairs(velInflow, Mxy_inflow/1000000, 'b-', 'LineWidth', 1.5);
ylabel('M_z [a.u.] (inflow enhancement)');
% ylim([0 1]);
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

if saveThis || ~exist(fullfile(info.project.figures,'radialProfiles.fig'),'file')
    saveas(        f, fullfile(info.project.figures,'radialProfiles.fig'));
    exportgraphics(f, fullfile(info.project.figures,'radialProfiles.png'));
    exportgraphics(f, fullfile(info.project.figures,'radialProfiles.svg'));
end
if saveThis || ~exist(fullfile(info.project.figures,'radialProfiles.mat'),'file')
    save(fullfile(info.project.figures,'radialProfiles.mat'), 'vel_fit', 'vel_fit2', 'mag_fit');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



if 0
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
I_ph    = squeeze(mean(data,[1 2]));
Ivenc_ph = squeeze(dataVenc);
Mnorm_ph = abs(mean(I_ph(Ivenc_ph==inf)));
vencList_ph = sort(unique(Ivenc_ph),'descend');
I_ph_mean = arrayfun(@(v) mean(I_ph(Ivenc_ph==v)), vencList_ph) / Mnorm_ph;
plotComplexDomain(ax{end}, I_ph_mean, 'full', 'markers');

% Save
if saveThis || ~exist(fullfile(info.project.figures,'phantomSummary.fig'))
    saveas(        f      , fullfile(info.project.figures,'phantomSummary.fig'));
    exportgraphics(f      , fullfile(info.project.figures,'phantomSummary.png'));
    exportgraphics(f      , fullfile(info.project.figures,'phantomSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% return

if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot matched simulation summary -- same format as phantom summary but all from simulation
%    spatial maps (rows=PE, cols=FE in runSim → transposed for display: rows=FE, cols=PE)
%    ISMRM2026-poster.pptx slide 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(info.project.figures,'radialProfiles.mat'), 'vel_fit', 'mag_fit');

% Set up simulation
p = runSim;
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

% match voxel grid to phantom data
pSim.voxGrid.fovFE = size(data,1) * FEspacing; % match phantom FOV exactly
pSim.voxGrid.fovPE = size(data,2) * PEspacing;
pSim.voxGrid.matFE = size(data,1);
pSim.voxGrid.matPE = size(data,2);
pSim.nSpin        = (2^10)^2;
pSim.gridMode     = 'pseudoVoxel';

% match vessel geometry to phantom data (fitted values)
pVessel.ID        = vel_fit.R*2;
pVessel.WT        = OD/2-ID/2; %2.38125;
% match vessel velocity to phantom data (fitted values again)
pVessel.vMean     = vel_fit.vMax / 2;
pVessel.profile   = 'parabolic1'; % with this, runSim.m should generate the velocity profile using the same function used for the velocity fit

% match MRI acquisition parameters
pMri.fieldStrength   = 3;
pMri.species         = 'phantom';
pMri.sliceThickness  = 2.2;              % [mm]
pMri.TR              = 75.90/(5+1)/1000; % [s]
pMri.TE              = 9.8/1000;         % [s]
pMri.FA              = 50;               % [deg]
pMri.venc.method     = 'FVEmono'; % use defaults

% precompute spinGrid using a dummy run
p = runSim(pVessel,pSim,pMri);
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

% match S to phantom data
[spinGridFE, spinGridPE] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);
spinGridR = sqrt(spinGridFE.^2 + spinGridPE.^2);
pVessel.S.lumen    = max(0, mag_fit(spinGridR(pVessel.mask.lumen))) ./pSim.nSpinPerVox;
pVessel.S.surround = mean(M(maskTissueOnly))                        ./pSim.nSpinPerVox;
pVessel.S.wall     = mean(M(maskWallOnly))                          ./pSim.nSpinPerVox;


% mag_fit

% figure
% imagesc(pVessel.mask.lumen)
% figure
% imagesc(spinGridR); colorbar
% figure
% imagesc(mag_fit(spinGridR))

% Run simulation
resSim = runSim(pVessel, pSim, pMri, [], false);

% Display variables
FEax_sim     = resSim.pSim.spinGrid.coorFE;   % [mm]
PEax_sim     = resSim.pSim.spinGrid.coorPE;   % [mm]
magMap_sim   = double(resSim.magMap);
vMap_sim     = double(resSim.vMap);
lumen_sim    = resSim.pVessel.mask.lumen;
wall_sim     = resSim.pVessel.mask.wall;
surround_sim = resSim.pVessel.mask.surround;
theta_sim    = linspace(0, 2*pi, 360);

% Complex-domain spiral: squeeze to [nVenc x 1], normalize by max magnitude
Ienc_sim      = squeeze(resSim.I);
Ienc_norm_sim = Ienc_sim ./ max(abs(Ienc_sim));

% Plot
fSim = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hTSim = tiledlayout(fSim,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); axSim = {};

% Magnitude map
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, magMap_sim, [0 max(magMap_sim(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end},'simulation ROI');

% Velocity map
axSim{end+1} = nexttile(hTSim);
vLim_sim = max(abs(vMap_sim(:)));
imagesc(axSim{end}, PEax_sim, FEax_sim, vMap_sim, [-vLim_sim vLim_sim]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axSim{end}.Colormap = redblue; set(axSim{end},'XTick',[],'YTick',[]);
title(axSim{end}, sprintf('parabolic vMean=%.1f cm/s', vel_fit.vMax/2));

% Lumen mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(lumen_sim), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on'); plot(axSim{end}, ID/2*cos(theta_sim), ID/2*sin(theta_sim), 'm');
title(axSim{end},'lumen mask');

% Wall mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(wall_sim), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on');
plot(axSim{end}, ID/2*cos(theta_sim),       ID/2*sin(theta_sim), 'm');
plot(axSim{end}, (ID/2+2.38125)*cos(theta_sim), (ID/2+2.38125)*sin(theta_sim), 'm');
title(axSim{end},'wall mask');

% Surround mask
axSim{end+1} = nexttile(hTSim);
imagesc(axSim{end}, PEax_sim, FEax_sim, single(surround_sim), [0 1]); axis image;
axSim{end}.Colormap = gray; set(axSim{end},'XTick',[],'YTick',[]);
hold(axSim{end},'on'); plot(axSim{end}, (ID/2+2.38125)*cos(theta_sim), (ID/2+2.38125)*sin(theta_sim), 'm');
title(axSim{end},'surround mask');

% Complex-plane spiral
axSim{end+1} = nexttile(hTSim, [3 3]);
plotComplexDomain(axSim{end}, Ienc_norm_sim, 'tight', 'line');
title(axSim{end}, 'complex-domain signal evolution (simulation)');

if saveThis || ~exist(fullfile(info.project.figures,'matchedSimSummary.fig'),'file')
    saveas(        fSim, fullfile(info.project.figures,'matchedSimSummary.fig'));
    exportgraphics(fSim, fullfile(info.project.figures,'matchedSimSummary.png'));
    exportgraphics(fSim, fullfile(info.project.figures,'matchedSimSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot matched simulation summary V2 -- poly4 mag profile + flat cap + vMean=4.648 cm/s
%    Based on multiVenc/doIt_phantom03_2.m (original slide 8 source)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(info.project.figures,'radialProfiles.mat'), 'vel_fit', 'mag_fit');

% Set up simulation
p = runSim;
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

% match voxel grid to phantom data
pSim.voxGrid.fovFE = size(data,1) * FEspacing;
pSim.voxGrid.fovPE = size(data,2) * PEspacing;
pSim.voxGrid.matFE = size(data,1);
pSim.voxGrid.matPE = size(data,2);
pSim.nSpin        = (2^10)^2;
pSim.gridMode     = 'pseudoVoxel';

% match vessel geometry — same as current
pVessel.ID      = vel_fit.R*2;
pVessel.WT      = OD/2-ID/2;
pVessel.profile = 'parabolic1';
% vMean from original doIt_phantom03_2.m (hardcoded, not from fit)
pVessel.vMean   = 4.648; % [cm/s]

% match MRI acquisition parameters
pMri.fieldStrength   = 3;
pMri.species         = 'phantom';
pMri.sliceThickness  = 2.2;
pMri.TR              = 75.90/(5+1)/1000;
pMri.TE              = 9.8/1000;
pMri.FA              = 50;
pMri.venc.method     = 'FVEmono';

% precompute spinGrid
p = runSim(pVessel, pSim, pMri);
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

% match S to phantom data — poly4 + flat cap (from doIt_phantom03_2.m lines 99-112)
mag_fit_v2 = fit([1:5]',[1:5]','poly4');
mag_fit_v2.p1 = -3.275935517412673e-08;
mag_fit_v2.p2 =  1.797301928434331e-07;
mag_fit_v2.p3 = -3.776838813872159e-07;
mag_fit_v2.p4 =  2.753845832568102e-07;
mag_fit_v2.p5 =  5.385039631953827e-07;
rPeak_v2 = 5.414525403728579e-01;  % [mm] flatten for r < rPeak

[spinGridFE2, spinGridPE2] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);
spinGridR2 = sqrt(spinGridFE2.^2 + spinGridPE2.^2);
r_lumen2   = spinGridR2(pVessel.mask.lumen);
S_lumen2   = mag_fit_v2(r_lumen2);
S_lumen2(r_lumen2 < rPeak_v2) = mag_fit_v2(rPeak_v2);  % flat cap
pVessel.S.lumen    = max(0, S_lumen2) ./ pSim.nSpinPerVox;
pVessel.S.surround = mean(M(maskTissueOnly)) ./ pSim.nSpinPerVox;
pVessel.S.wall     = mean(M(maskWallOnly))   ./ pSim.nSpinPerVox;

% Run simulation
resSim2 = runSim(pVessel, pSim, pMri, [], false);

% Display variables
FEax_sim2     = resSim2.pSim.spinGrid.coorFE;
PEax_sim2     = resSim2.pSim.spinGrid.coorPE;
magMap_sim2   = double(resSim2.magMap);
vMap_sim2     = double(resSim2.vMap);
lumen_sim2    = resSim2.pVessel.mask.lumen;
wall_sim2     = resSim2.pVessel.mask.wall;
surround_sim2 = resSim2.pVessel.mask.surround;
theta_sim2    = linspace(0, 2*pi, 360);
Ienc_sim2      = squeeze(resSim2.I);
Ienc_norm_sim2 = Ienc_sim2 ./ max(abs(Ienc_sim2));

% Plot
fSim2 = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hTSim2 = tiledlayout(fSim2,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); axSim2 = {};

axSim2{end+1} = nexttile(hTSim2);
imagesc(axSim2{end}, PEax_sim2, FEax_sim2, magMap_sim2, [0 max(magMap_sim2(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axSim2{end}.Colormap = gray; set(axSim2{end},'XTick',[],'YTick',[]);
title(axSim2{end},'simulation ROI');

axSim2{end+1} = nexttile(hTSim2);
vLim_sim2 = max(abs(vMap_sim2(:)));
imagesc(axSim2{end}, PEax_sim2, FEax_sim2, vMap_sim2, [-vLim_sim2 vLim_sim2]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axSim2{end}.Colormap = redblue; set(axSim2{end},'XTick',[],'YTick',[]);
title(axSim2{end}, 'parabolic vMean=4.648 cm/s');

axSim2{end+1} = nexttile(hTSim2);
imagesc(axSim2{end}, PEax_sim2, FEax_sim2, single(lumen_sim2), [0 1]); axis image;
axSim2{end}.Colormap = gray; set(axSim2{end},'XTick',[],'YTick',[]);
hold(axSim2{end},'on'); plot(axSim2{end}, ID/2*cos(theta_sim2), ID/2*sin(theta_sim2), 'm');
title(axSim2{end},'lumen mask');

axSim2{end+1} = nexttile(hTSim2);
imagesc(axSim2{end}, PEax_sim2, FEax_sim2, single(wall_sim2), [0 1]); axis image;
axSim2{end}.Colormap = gray; set(axSim2{end},'XTick',[],'YTick',[]);
hold(axSim2{end},'on');
plot(axSim2{end}, ID/2*cos(theta_sim2), ID/2*sin(theta_sim2), 'm');
plot(axSim2{end}, (ID/2+2.38125)*cos(theta_sim2), (ID/2+2.38125)*sin(theta_sim2), 'm');
title(axSim2{end},'wall mask');

axSim2{end+1} = nexttile(hTSim2);
imagesc(axSim2{end}, PEax_sim2, FEax_sim2, single(surround_sim2), [0 1]); axis image;
axSim2{end}.Colormap = gray; set(axSim2{end},'XTick',[],'YTick',[]);
hold(axSim2{end},'on'); plot(axSim2{end}, (ID/2+2.38125)*cos(theta_sim2), (ID/2+2.38125)*sin(theta_sim2), 'm');
title(axSim2{end},'surround mask');

axSim2{end+1} = nexttile(hTSim2, [3 3]);
plotComplexDomain(axSim2{end}, Ienc_norm_sim2, 'tight', 'line');
title(axSim2{end}, 'complex-domain signal (poly4+cap, vMean=4.648)');

if saveThis || ~exist(fullfile(info.project.figures,'matchedSimSummaryV2.fig'),'file')
    saveas(        fSim2, fullfile(info.project.figures,'matchedSimSummaryV2.fig'));
    exportgraphics(fSim2, fullfile(info.project.figures,'matchedSimSummaryV2.png'));
    exportgraphics(fSim2, fullfile(info.project.figures,'matchedSimSummaryV2.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


if 0
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


if 0
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot in vivo summary -- sub-01 vessel-01, reference mag, velocity map, complex-domain signal
%    on the model of phantom summary -- ISMRM2026-poster.pptx slide 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 1; roiIdx = 1;
img      = inVivoSubData{s}.img;
imgInfo  = inVivoSubData{s}.imgInfo;
refImgAv = inVivoSubData{s}.refImgAv;
roiY     = inVivoSubRoiList{s}(roiIdx).roiY;
roiX     = inVivoSubRoiList{s}(roiIdx).roiX;
subName  = inVivoSubNames{s};

% Reference magnitude (full slice)
refMag_iv = squeeze(mean(abs(refImgAv), 7));

% Velocity map at middle venc (averaged over runs and repetitions)
vencList_iv = sort(unique(imgInfo.vencList(isfinite(imgInfo.vencList))));
vencSel_iv  = vencList_iv(round(end/2));
refAvg_iv   = squeeze(mean(img(:,:,:,:,:,:,imgInfo.vencList==inf,         :,:,:,:,:,:,:,:,:), [7 11]));
selAvg_iv   = squeeze(mean(img(:,:,:,:,:,:,imgInfo.vencList==vencSel_iv,  :,:,:,:,:,:,:,:,:), [7 11]));
velMap_iv   = phase2vel(angle(selAvg_iv ./ refAvg_iv), vencToM1(vencSel_iv));

% ROI overlay -- separate figure
hFroi = figure('MenuBar','none','ToolBar','none','Units','normalized','Position',[0 0 1 1]);
hIm_roi = imagesc(refMag_iv);
hold on;
plot([roiX(1)-.5 roiX(2)+.5 roiX(2)+.5 roiX(1)-.5 roiX(1)-.5], ...
     [roiY(1)-.5 roiY(1)-.5 roiY(2)+.5 roiY(2)+.5 roiY(1)-.5], 'c', 'LineWidth', 1);
text(roiX(1), roiY(1), 'vessel01', 'Color','r','FontSize',12,'FontWeight','bold');
axRoi = gca; axis image;
set(axRoi,'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
axRoi.XAxis.Visible = 'off'; axRoi.YAxis.Visible = 'off';
axRoi.CLim = [0, 0.5*max(hIm_roi.CData(:))];
figName_roi = [subName '-vessel01-roiOverlay'];
if saveThis || ~exist(fullfile(info.project.figures,[figName_roi '.png']),'file')
    drawnow;
    exportgraphics(hFroi, fullfile(info.project.figures,[figName_roi '.png']));
    exportgraphics(hFroi, fullfile(info.project.figures,[figName_roi '.svg']));
end
close(hFroi);

% Extract complex-domain signal for this ROI
trjIV = img(roiY(1):roiY(2), roiX(1):roiX(2), :,:,:,:,:,:,:,:,:,:,:,:,:,:);
runIdxList_iv = unique(imgInfo.runIdx);
for rr = 1:length(runIdxList_iv)
    idx_rr  = squeeze(imgInfo.runIdx==runIdxList_iv(rr) & imgInfo.vencList==inf);
    refPhaseIV = angle(mean(trjIV(:,:,:,:,:,:,idx_rr,:,:,:,:,:,:,:,:,:), [7 11]));
    idx_rr2 = squeeze(imgInfo.runIdx==runIdxList_iv(rr));
    trjIV(:,:,:,:,:,:,idx_rr2,:,:,:,:,:,:,:,:,:) = ...
        trjIV(:,:,:,:,:,:,idx_rr2,:,:,:,:,:,:,:,:,:) ./ exp(1i*refPhaseIV);
end
trjIV      = permute(mean(trjIV,[1 2]),[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
trjIV      = trjIV ./ abs(mean(trjIV(1:2,:),[1 2]));
trjVencIV  = permute(imgInfo.vencList,[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);

% Summary figure -- 3 rows × 4 cols, column-major
f_ivs = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT_ivs = tiledlayout(f_ivs,3,4,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor');
ax_ivs = {};

% Full-slice reference magnitude with ROI box
ax_ivs{end+1} = nexttile(hT_ivs);
imagesc(ax_ivs{end}, refMag_iv, [0, 0.5*max(refMag_iv(:))]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
hold(ax_ivs{end},'on');
plot(ax_ivs{end}, [roiX(1)-.5 roiX(2)+.5 roiX(2)+.5 roiX(1)-.5 roiX(1)-.5], ...
                  [roiY(1)-.5 roiY(1)-.5 roiY(2)+.5 roiY(2)+.5 roiY(1)-.5], 'c', 'LineWidth', 1);
title(ax_ivs{end}, [subName ' vessel 01']);

% ROI crop -- reference magnitude
ax_ivs{end+1} = nexttile(hT_ivs);
refMagCrop_iv = refMag_iv(roiY(1):roiY(2), roiX(1):roiX(2));
imagesc(ax_ivs{end}, refMagCrop_iv, [0, max(refMagCrop_iv(:))]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
ylabel(colorbar(ax_ivs{end},'Location','westoutside'), 'MR magn. [a.u.]');
title(ax_ivs{end}, 'ROI mag');

% ROI crop -- velocity map
ax_ivs{end+1} = nexttile(hT_ivs);
velCrop_iv = velMap_iv(roiY(1):roiY(2), roiX(1):roiX(2));
vLim_ivs = max(abs(velCrop_iv(:)));
imagesc(ax_ivs{end}, velCrop_iv, [-vLim_ivs vLim_ivs]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',redblue,'DataAspectRatio',[imgInfo.res 1]);
ylabel(colorbar(ax_ivs{end},'Location','westoutside'), 'velocity [cm/s]');
title(ax_ivs{end}, sprintf('venc=%gcm/s', vencSel_iv));

% Complex-domain signal (3 rows × 3 cols span)
ax_ivs{end+1} = nexttile(hT_ivs,[3 3]);
plotComplexDomain(ax_ivs{end}, trjIV(:), 'tight', 'markers');

% Save
figName_ivs = [subName '-vessel01-summary'];
if saveThis || ~exist(fullfile(info.project.figures,[figName_ivs '.fig']),'file')
    saveas(        f_ivs, fullfile(info.project.figures,[figName_ivs '.fig']));
    exportgraphics(f_ivs, fullfile(info.project.figures,[figName_ivs '.png']));
    exportgraphics(f_ivs, fullfile(info.project.figures,[figName_ivs '.svg']));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

return

if 0
saveThis = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot in vivo -- reference mag, velocity map and complex-domain signal evolution -- ISMRM2026-poster.pptx slide 10 (sub-01) and 11 (sub-02)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(info.project.figures,'dir'); mkdir(info.project.figures); end

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
    if saveThis || ~exist(fullfile(info.project.figures,[inVivoSubNames{s} '-roiOverlay.png']),'file')
        drawnow;
        exportgraphics(hF, fullfile(info.project.figures,[inVivoSubNames{s} '-roiOverlay.png']));
        exportgraphics(hF, fullfile(info.project.figures,[inVivoSubNames{s} '-roiOverlay.svg']));
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
        if saveThis || ~exist(fullfile(info.project.figures,[vesselFigName '.png']),'file')
            drawnow;
            exportgraphics(hFv, fullfile(info.project.figures,[vesselFigName '.png']));
            exportgraphics(hFv, fullfile(info.project.figures,[vesselFigName '.svg']));
        end
        close(hFv); clear hPcont d prob
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




if 0
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
[N,edges] = histcounts(res.vMap(getVoxIdx(res.pSim.voxGrid,res.pSim.spinGrid)==0),20);
% binWidth = mean(diff(edges));
% edges = edges-binWidth/2; edges(end+1) = edges(end)+binWidth;
% [N,edges] = histcounts(res.vMap(getVoxIdx(res.pSim.voxGrid,res.pSim.spinGrid)==0),edges);
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

if ~exist(info.project.figures,'dir'); mkdir(info.project.figures); end
if saveThis || ~exist(fullfile(info.project.figures,'FVEvelSpec.fig'),'file') || ~exist(fullfile(info.project.figures,'FVEvelSpec_inflow.fig'),'file')
    saveas(        fVelSpec      , fullfile(info.project.figures,'FVEvelSpec.fig'       ));
    exportgraphics(fVelSpec      , fullfile(info.project.figures,'FVEvelSpec.png'       ));
    exportgraphics(fVelSpec      , fullfile(info.project.figures,'FVEvelSpec.svg'       ));
    saveas(        fVelSpecInflow, fullfile(info.project.figures,'FVEvelSpec_inflow.fig'));
    exportgraphics(fVelSpecInflow, fullfile(info.project.figures,'FVEvelSpec_inflow.png'));
    exportgraphics(fVelSpecInflow, fullfile(info.project.figures,'FVEvelSpec_inflow.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVelSpec; % FVE spectra reflects spin velocity distribution, but weighted by velocity-dependent spin magnitude
fVelSpecInflow; % Here the weighting effect was maximized using a 90 flip angle for a linear magnitude function of velocity
end







