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

%Scale magnitude to a decent range
magScaleFac = mean(abs(data(:)));
data       = data      ./magScaleFac;
dataNoFlow = dataNoFlow./magScaleFac;


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
%%%%%%%%%%%%%%%%%%%%%%%
%% Radial profiles fits
%%%%%%%%%%%%%%%%%%%%%%%
cFlow          = mean(data(:,:,dataVenc==inf),3);
cFlowSpoiled   = mean(data(:,:,dataVenc==min(dataVenc)),3);
cNoFlow        = mean(dataNoFlow(:,:,dataVenc==inf),3);
cNoFlowSpoiled = mean(dataNoFlow(:,:,dataVenc==min(dataVenc)),3);
vFlow   = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
vNoFlow = zeros(size(vFlow));
rFlow   = rGrid;
rNoFlow = rGrid;

% Fit of v(r) velocity function of radius
R = ID/2;
[velFit, velFitFixR] = fitVelProfile(rGrid(maskBloodOnly), vFlow(maskBloodOnly), R);

% Fit of m(r) magnitude function of radius
R = ID/2;
B = abs(mean(cNoFlow(maskBloodOnly)));
[magFit, magFitFixR, magFitFixB, magFitFixBR] = fitMagProfile(rGrid(maskBloodOnly), abs(cFlow(maskBloodOnly)), R, B);

% Fit of m(r) magnitude function of radius for the spoiled low-venc data
R = ID/2;
B = abs(mean(cNoFlowSpoiled(maskBloodOnly)));
[magSpoilFit, magSpoilFitFixR, magSpoilFitFixB, magSpoilFitFixBR] = fitMagProfile(rGrid(maskBloodOnly), abs(cFlowSpoiled(maskBloodOnly)), R, B);

% Fit of v(r) and m(v(r)) jointly
R = ID/2;
[velJointFit, magJoinFit] = fitMagVelProfile(rGrid(maskBloodOnly), vFlow(maskBloodOnly), abs(cFlow(maskBloodOnly)), abs(cNoFlow(maskBloodOnly)), R, 'joint', 4);


f = figure;
hT = tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact'); ax = {};
theta = linspace(0, 2*pi, 360);

%mag
ax{end+1} = nexttile;
imagesc(PEpos,FEpos,abs(cFlow),[0 max(abs(cFlow(maskBloodOnly)))]); axis image; colormap(ax{end},'gray'); colorbar;
hold on
plot(R*cos(theta), R*sin(theta), 'w--', 'LineWidth', 1);

ax{end+1} = nexttile;
imagesc(PEpos,FEpos,reshape(magJoinFit(velJointFit(rGrid)),size(rGrid))); axis image; colormap(ax{end},'gray'); colorbar;
hold on
plot(R*cos(theta), R*sin(theta), 'w--', 'LineWidth', 1);

cLim = get([ax{end-1:end}],'CLim'); cLim = [0 1] .*max([cLim{:}]); set([ax{end-1:end}],'CLim',cLim);


%vel
ax{end+1} = nexttile;
imagesc(PEpos,FEpos,vFlow,[-1 1].*max(vFlow(maskBloodOnly))); axis image; colormap(ax{end},'jet'); colorbar;
hold on
plot(R*cos(theta), R*sin(theta), 'w--', 'LineWidth', 1);

ax{end+1} = nexttile;
imagesc(PEpos,FEpos,reshape(velJointFit(rGrid),size(rGrid))); axis image; colormap(ax{end},'jet'); colorbar;
hold on
plot(R*cos(theta), R*sin(theta), 'w--', 'LineWidth', 1);

cLim = get([ax{end-1:end}],'CLim'); cLim = [-1 1] .*max(abs([cLim{:}])); set([ax{end-1:end}],'CLim',cLim);


%mag vs vel
ax{end+1} = nexttile;
scatter(vFlow(maskBloodOnly),abs(cFlow(maskBloodOnly)),'MarkerFaceColor','w','MarkerEdgeColor','k')
axis square; xlabel('velocity'); ylabel('magnitude'); hold on
v = linspace(0,max(vFlow(maskBloodOnly)),100);
plot(v,magJoinFit(v),'w-');
grid on

ax{end+1} = nexttile;
scatter(velJointFit(rGrid(maskBloodOnly)),magJoinFit(velJointFit(rGrid(maskBloodOnly))),'MarkerFaceColor','w','MarkerEdgeColor','k')
axis square; xlabel('velocity'); ylabel('magnitude'); hold on
plot(v,magJoinFit(v),'w-');
grid on

xLim = get([ax{end-1:end}],'XLim'); xLim = [0 1] .*max([xLim{:}]); set([ax{end-1:end}],'XLim',xLim);
yLim = get([ax{end-1:end}],'YLim'); yLim = [0 1] .*max([yLim{:}]); set([ax{end-1:end}],'YLim',yLim);

if saveThis || ~exist(fullfile(info.project.figures,'radialProfilesFits.fig'),'file')
    saveas(        f, fullfile(info.project.figures,'radialProfilesFits.fig'));
    exportgraphics(f, fullfile(info.project.figures,'radialProfilesFits.png'));
    exportgraphics(f, fullfile(info.project.figures,'radialProfilesFits.svg'));
end
if saveThis || ~exist(fullfile(info.project.figures,'radialProfilesFits.mat'),'file')
    save(fullfile(info.project.figures,'radialProfilesFits.mat'), ...
        'velFit', 'velFitFixR', ...
        'magFit', 'magFitFixR', 'magFitFixB', 'magFitFixBR', ...
        'magSpoilFit', 'magSpoilFitFixR', 'magSpoilFitFixB', 'magSpoilFitFixBR', ...
        'velJointFit', 'magJoinFit');
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



if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot matched simulation summary — joint velocity + magnitude fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(info.project.figures,'radialProfilesFits.mat'), 'velJointFit', 'magJoinFit');

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
pSim.nSpin         = (2^10)^2;
pSim.gridMode      = 'pseudoVoxel';

% match vessel geometry from joint velocity fit
pVessel.ID      = velJointFit.R * 2;
pVessel.WT      = OD/2 - ID/2;
pVessel.vMean   = velJointFit.Vmax / 2;
pVessel.profile = 'parabolic1';

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

% match S to phantom data using joint m(v(r)) fit
[spinGridFE3, spinGridPE3] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);
spinGridR3 = sqrt(spinGridFE3.^2 + spinGridPE3.^2);
r_lumen3   = spinGridR3(pVessel.mask.lumen);
pVessel.S.lumen    = max(0, magJoinFit(velJointFit(r_lumen3)))    ./ pSim.nSpinPerVox;
pVessel.S.surround = mean(abs(cFlow(maskTissueOnly)))              ./ pSim.nSpinPerVox;
pVessel.S.wall     = mean(abs(cFlow(maskWallOnly)))                ./ pSim.nSpinPerVox;

% Run simulation
resSim3 = runSim(pVessel, pSim, pMri, [], false);

% Display variables
FEax3      = resSim3.pSim.spinGrid.coorFE;
PEax3      = resSim3.pSim.spinGrid.coorPE;
magMap3    = double(resSim3.magMap);
vMap3      = double(resSim3.vMap);
theta3     = linspace(0, 2*pi, 360);
Ienc3      = squeeze(resSim3.I);
Ienc3_norm = Ienc3 ./ max(abs(Ienc3));

% Plot
fSim3  = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hTSim3 = tiledlayout(fSim3,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); axSim3 = {};

axSim3{end+1} = nexttile(hTSim3);
imagesc(axSim3{end}, PEax3, FEax3, magMap3, [0 max(magMap3(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axSim3{end}.Colormap = gray; set(axSim3{end},'XTick',[],'YTick',[]);
title(axSim3{end}, 'simulation ROI');

axSim3{end+1} = nexttile(hTSim3);
imagesc(axSim3{end}, PEax3, FEax3, vMap3, [-1 1].*max(abs(vMap3(:)))); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axSim3{end}.Colormap = redblue; set(axSim3{end},'XTick',[],'YTick',[]);
title(axSim3{end}, sprintf('joint fit: vMean=%.2f cm/s, R=%.2f mm', velJointFit.Vmax/2, velJointFit.R));

axSim3{end+1} = nexttile(hTSim3);
imagesc(axSim3{end}, PEax3, FEax3, single(resSim3.pVessel.mask.lumen), [0 1]); axis image;
axSim3{end}.Colormap = gray; set(axSim3{end},'XTick',[],'YTick',[]);
hold(axSim3{end},'on'); plot(axSim3{end}, ID/2*cos(theta3), ID/2*sin(theta3), 'm');
title(axSim3{end}, 'lumen mask');

axSim3{end+1} = nexttile(hTSim3);
imagesc(axSim3{end}, PEax3, FEax3, single(resSim3.pVessel.mask.wall), [0 1]); axis image;
axSim3{end}.Colormap = gray; set(axSim3{end},'XTick',[],'YTick',[]);
hold(axSim3{end},'on');
plot(axSim3{end}, ID/2*cos(theta3),         ID/2*sin(theta3),         'm');
plot(axSim3{end}, (OD/2)*cos(theta3), (OD/2)*sin(theta3), 'm');
title(axSim3{end}, 'wall mask');

axSim3{end+1} = nexttile(hTSim3);
imagesc(axSim3{end}, PEax3, FEax3, single(resSim3.pVessel.mask.surround), [0 1]); axis image;
axSim3{end}.Colormap = gray; set(axSim3{end},'XTick',[],'YTick',[]);
hold(axSim3{end},'on'); plot(axSim3{end}, (OD/2)*cos(theta3), (OD/2)*sin(theta3), 'm');
title(axSim3{end}, 'surround mask');

axSim3{end+1} = nexttile(hTSim3, [3 3]);
plotComplexDomain(axSim3{end}, Ienc3_norm, 'tight', 'line');
title(axSim3{end}, 'complex-domain signal (joint v(r)+m(v) fit)');

if saveThis || ~exist(fullfile(info.project.figures,'matchedSimSummaryFits.fig'),'file')
    saveas(        fSim3, fullfile(info.project.figures,'matchedSimSummaryFits.fig'));
    exportgraphics(fSim3, fullfile(info.project.figures,'matchedSimSummaryFits.png'));
    exportgraphics(fSim3, fullfile(info.project.figures,'matchedSimSummaryFits.svg'));
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







