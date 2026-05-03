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
    dbstack; error('double check that');
    envId = 1;
    storageDrive  = '/local/users/Proulx-S/';
    scratchDrive  = '/scratch/users/Proulx-S/';
    databaseDrive = '/local/users/Proulx-S/db/';
    databasePhantomDrive = '/local/users/Proulx-S/dbPhantom/';
    projectCode    = fullfile(scratchDrive, projectName);        if ~exist(projectCode,'dir');    mkdir(projectCode);    end
    projectStorage = fullfile(storageDrive, projectName);        if ~exist(projectStorage,'dir'); mkdir(projectStorage); end
    projectScratch = fullfile(scratchDrive, projectName, 'tmp'); if ~exist(projectScratch,'dir'); mkdir(projectScratch); end
    toolDir        = '/scratch/users/Proulx-S/tools';            if ~exist(toolDir,'dir');        mkdir(toolDir);        end
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
    projectDataBase = fullfile(databaseDrive                               ); if ~exist(projectDataBase,'dir'); mkdir(projectDataBase); end
    projectDataBasePhantom = fullfile(databasePhantomDrive          ); if ~exist(projectDataBasePhantom,'dir'); mkdir(projectDataBasePhantom); end
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
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Illustrate the effect of inflow on FVE spectrum
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



forceThis = 0;
%%%%%%%%%%%%%%%%%%%%
%% Load phantom data
%%%%%%%%%%%%%%%%%%%%
phantom03dataFile = fullfile(info.project.scratch, 'phantom03.mat');
if forceThis || ~exist(phantom03dataFile,'file')
    [data, dataVenc, dataRun, dataNoFlow, PEspacing, FEspacing] = loadPhantom03(fullfile(info.project.dataBasePhantom,'20251010_multiVENCphantom03'));
    save(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
else
    load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
end

ID =  6.35; % mm
OD = 11.11; % mm

data = conj(data);
dataNoFlow = conj(dataNoFlow);

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
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom summary -- reference mag, velocity map, masks, and complex-domain signal evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT = tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};

% Plot mag
ax{end+1} = nexttile;
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
hIm = imagesc(ax{end},PEpos,FEpos,M,[0 max(M(:))]); axis image;
% hIm = imagesc(M,[0 max(M(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray;
set(ax{end},'XTick',[],'YTick',[]);
% set(ax{end},'DataAspectRatio',[FEspacing/PEspacing 1 1],'XTick',[],'YTick',[]);
title(ax{end},'phantom ROI');

% Plot velocity map
ax{end+1} = nexttile;
dataVencIdx    = find(dataVenc==10);
dataVencRefIdx = find(dataVenc==inf);
venc = dataVenc(dataVencIdx); venc = unique(venc);
PD   = angle(mean(  data(:,:,dataVencIdx) ./ exp(1j.*angle(mean(data(:,:,dataVencRefIdx),3)))  ,3));
CD   = mean(data(:,:,dataVencIdx),3)-mean(data(:,:,dataVencRefIdx),3);
[velCD,phi,velPD,~,~,~] = getPlugFlowEstimates(venc,CD,[],[],PD);
hIm = imagesc(ax{end},PEpos,FEpos,-velPD,[-max(abs(velPD(:))) max(abs(velPD(:)))]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = redblue;
set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},['venc=' num2str(venc) 'cm/s']);

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
theta = linspace(0, 2*pi, 360);                                                                                                               
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

if 0
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot phantom details -- all maps with and without flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mag flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};
vencList = sort(unique(dataVenc),'descend'); M = {};
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
    P{end+1} = -angle(mean(data(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,P{end},[-pi pi]); axis image
    title(ax{end},['flow on; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'MR magn. [a.u.]');
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
    P{end+1} = -angle(mean(dataNoFlow(:,:,dataVenc==vencList(vencIdx)),3));
    imagesc(ax{end},PEpos,FEpos,P{end},[-pi pi]); axis image
    title(ax{end},['flow off; venc=' num2str(vencList(vencIdx)) ' cm/s']);
    set(ax{end},'XTick',[],'YTick',[]);
    if vencIdx<5
        cb = colorbar(ax{end},'Location','westoutside');
        ylabel(cb, 'MR magn. [a.u.]');
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

return


% p = runSim;
% tmp = runSim(p.pVessel,p.pSim,p.pMri);
% p.pVessel.vMean = tmp.pMri.vCrit/2*1.5;
% p.pMri.FA = 90;
% p.pMri.venc.method = 'FVEmono';
% p.pMri.venc.FVEbw = p.pVessel.vMean*4;
% p.pMri.venc.FVEres = p.pMri.venc.FVEbw./200;
% p.pSim.nSpin = (2^10+1)^2;
% res      = runSim(p.pVessel,p.pSim,p.pMri,[],0);
% % p.pMri.sliceThickness = inf;
% Mz  = getMz_ss(p.pMri,p.pMri.relax.blood,p.pVessel.vMean);
% Mxy = getMxy_ss(Mz,p.pMri,p.pMri.relax.blood);
% p.pVessel.S.lumen = Mxy;
% resSatin = runSim(p.pVessel,p.pSim,p.pMri);



% figure('MenuBar','none','ToolBar','none');
% imagesc(tmp1>0.95e-7,[0 max(tmp1(:))]); axis image;
% scatter(PD(tmp1>0.95e-7),phi(tmp1>0.95e-7)); axis image; grid on;
% xlabel('PD'); ylabel('phi');
% lim = [xlim ylim]; lim = [min(lim) max(lim)];
% hold on
% line(lim,lim)

% figure('MenuBar','none','ToolBar','none');
% scatter(PD(tmp1>0.95e-7),tmp1(tmp1>0.95e-7)); grid on;
% scatter(PD(:),tmp1(:)); grid on;
% xlabel('PD'); ylabel('mag');
% lim = [xlim ylim]; lim = [min(lim) max(lim)];
% hold on
% line(lim,lim)










% plug flow simulation
pVesselPlug = pVessel;
pVesselPlug.PD = pVesselPlug.ID/1;
resPlug = runSim(pVesselPlug, pVenc, pSim, pReal);

resTmp = resPlug;
ax{end+1} = nexttile;
tmpMagMap = resTmp.magMap;
imagesc(tmpMagMap); axis image;
ax{end}.Colormap = gray;
ax{end}.CLim = [0 max(tmpMagMap(:))*1.2];
set(ax{end},'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
title(ax{end},{'simulated ROI','plug flow'});

ax{end+1} = nexttile;
imagesc(resTmp.vMap); axis image;
ax{end}.Colormap = jet;
ax{end}.CLim = ax{2}.CLim;
set(ax{end},'XTick',[],'YTick',[]);

% laminar flow simulation
resTmp = res;
ax{end+1} = nexttile;
imagesc(tmpMagMap); axis image;
ax{end}.Colormap = gray;
ax{end}.CLim = [0 max(tmpMagMap(:))*1.2];
set(ax{end},'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
title(ax{end},{'simulated ROI','parabolic laminar flow'});

ax{end+1} = nexttile;
imagesc(resTmp.vMap); axis image;
ax{end}.Colormap = jet;
ax{end}.CLim = ax{2}.CLim;
set(ax{end},'XTick',[],'YTick',[]);

% blunted laminar flow simulation
pVesselBlunted = pVessel;
pVesselBlunted.PD = pVesselBlunted.ID/2;
resBlunted = runSim(pVesselBlunted, pVenc, pSim, pReal);

resTmp = resBlunted;
ax{end+1} = nexttile;
imagesc(tmpMagMap); axis image;
ax{end}.Colormap = gray;
ax{end}.CLim = [0 max(tmpMagMap(:))*1.2];
set(ax{end},'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
title(ax{end},{'simulated ROI','blunted laminar flow'});

ax{end+1} = nexttile;
imagesc(resTmp.vMap); axis image;
ax{end}.Colormap = jet;
ax{end}.CLim = ax{2}.CLim;
set(ax{end},'XTick',[],'YTick',[]);
scaleFont(f,1.2);

saveas(f,fullfile(pwd,'phantom03_phantomAndSimulationSetup.png'))
saveas(f,fullfile(pwd,'phantom03_phantomAndSimulationSetup.fig'))
saveas(f,fullfile(pwd,'phantom03_phantomAndSimulationSetup.eps'))















