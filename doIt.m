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


return


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
[inflowMz,~,~,~,~,inflowVel] = getMz_ss(p.pMri,p.pMri.relax.blood);
% inflowMxy = getMxy_ss(inflowMz,p.pMri,p.pMri.relax.blood);
hStairs = stairs(inflowVel,inflowMz,'g');
axis tight square; grid on; xlabel('spin velocity (cm/s)'); ylabel('M_z');
ylim([0 1])

if ~exist(fullfile(projectStorage, 'figures'),'dir'); mkdir(fullfile(projectStorage, 'figures')); end
if 0
    saveas(        fVelSpec      , fullfile(projectStorage, 'figures', 'FVEvelSpec.fig'       )                   );
    exportgraphics(fVelSpec      , fullfile(projectStorage, 'figures', 'FVEvelSpec.png'       ), 'Resolution', 300);
    exportgraphics(fVelSpec      , fullfile(projectStorage, 'figures', 'FVEvelSpec.svg'       )                   );
    saveas(        fVelSpecInflow, fullfile(projectStorage, 'figures', 'FVEvelSpec_inflow.fig')                   );
    exportgraphics(fVelSpecInflow, fullfile(projectStorage, 'figures', 'FVEvelSpec_inflow.png'), 'Resolution', 300);
    exportgraphics(fVelSpecInflow, fullfile(projectStorage, 'figures', 'FVEvelSpec_inflow.svg')                   );
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fVelSpec; % FVE spectra reflects spin velocity distribution, but weighted by velocity-dependent spin magnitude
fVelSpecInflow; % Here the weighting effect was maximized using a 90 flip angle for a linear magnitude function of velocity



%%%%%%%%%%%%%%%%%%%%
%% Load phantom data
%%%%%%%%%%%%%%%%%%%%
phantom03dataFile = fullfile(info.project.scratch, 'phantom03.mat');
if ~exist(phantom03dataFile,'file')
    [data, dataVenc,PEspacing,FEspacing] = loadPhantom03(fullfile(info.project.dataBasePhantom,'20251010_multiVENCphantom03'));
    save(phantom03dataFile, 'data', 'dataVenc', 'PEspacing', 'FEspacing');
else
    load(phantom03dataFile, 'data', 'dataVenc', 'PEspacing', 'FEspacing');
end
data = data ./ exp(1j*angle(data(dataVenc==inf,:,:)));
%% %%%%%%%%%%%%%%%%%




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



figure('MenuBar','none','ToolBar','none');
imagesc(tmp1>0.95e-7,[0 max(tmp1(:))]); axis image;
scatter(PD(tmp1>0.95e-7),phi(tmp1>0.95e-7)); axis image; grid on;
xlabel('PD'); ylabel('phi');
lim = [xlim ylim]; lim = [min(lim) max(lim)];
hold on
line(lim,lim)

figure('MenuBar','none','ToolBar','none');
scatter(PD(tmp1>0.95e-7),tmp1(tmp1>0.95e-7)); grid on;
scatter(PD(:),tmp1(:)); grid on;
xlabel('PD'); ylabel('mag');
lim = [xlim ylim]; lim = [min(lim) max(lim)];
hold on
line(lim,lim)




% Fig1A
f = figure('MenuBar','none','ToolBar','none','Units','Centimeter','Position',[7     9    24    10]);
hT = tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};

% phantom data
ax{end+1} = nexttile;
dataVencIdx = dataVenc==inf;
tmp1 = squeeze(abs(data(dataVencIdx,:,:)));
imagesc(tmp1,[0 max(tmp1(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray;
set(ax{end},'DataAspectRatio',[FEspacing/PEspacing 1 1],'XTick',[],'YTick',[]);
title(ax{end},'phantom ROI');

ax{end+1} = nexttile;
dataVencIdx    = find(dataVenc==11);
dataVencRefIdx = find(dataVenc==inf);
venc = dataVenc(dataVencIdx);
PD   = angle(data(dataVencIdx,:,:));
CD   = data(dataVencIdx,:,:)-data(dataVencRefIdx,:,:);
[velCD,phi,velPD,~,~,~] = getPlugFlowEstimates(venc,CD,[],[],PD);
imagesc(squeeze(velPD),[-max(velPD(:)) max(velPD(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = jet;
set(ax{end},'DataAspectRatio',[FEspacing/PEspacing 1 1],'XTick',[],'YTick',[]);



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















