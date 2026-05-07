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
% blueBlackRed (in util) auto-downloads Colorspace-Transformations on first call
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




if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 - motivation and goal
%%%%%%%%%%%%%%%%%%%%%%%%%%
% No MATLAB output — poster section only.
%% %%%%%%%%%%%%%%%%%%%%
end




if 1
%%%%%%%%%%%%%%%%%%%%%%%
%% 2 - hamilton context
%%%%%%%%%%%%%%%%%%%%%%%
% No MATLAB output — poster section only.
%% %%%%%%%%%%%%%%%%%%%%%
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

% %Scale magnitude to a decent range
% magScaleFac = mean(abs(data(:)));
% data       = data      ./magScaleFac;
% dataNoFlow = dataNoFlow./magScaleFac;


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
pGrid = -atan2(FEgrid, PEgrid);  % theta=0 → +PE ("right" in imagesc display)

theta = linspace(0, 2*pi, 360);

clear dFE dPE d_far d_near M com total
%% %%%%%%%%%%%%%%%%%



if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3 - education simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec3fig = fullfile(info.project.figures, '3-education-simulations');
if ~exist(sec3fig,'dir'); mkdir(sec3fig); end

pEdu = runSim; % get default parameters
pEdu.pSim.voxGrid.fovFE = size(data,1) * FEspacing;
pEdu.pSim.voxGrid.fovPE = size(data,2) * PEspacing;
pEdu.pSim.voxGrid.matFE = size(data,1);
pEdu.pSim.voxGrid.matPE = size(data,2);
pEdu.pSim.nSpin         = (2^10)^2;
pEdu.pSim.gridMode      = 'pseudoVoxel';
pEdu.pVessel.ID    = ID;
pEdu.pVessel.vMean = 5;
pEdu.pMri.fieldStrength = 3;
pEdu.pMri.species       = 'phantom';
pEdu.pMri.venc.FVEbw = 200;

pEdu = runSim(pEdu.pVessel,pEdu.pSim,pEdu.pMri); % update paremeters


% Flat magnitude profile: S.lumen = constant Mxy at vMean (no velocity-dependent inflow saturation)
Mz_flat  = getMz_ss(pEdu.pMri, pEdu.pMri.relax.blood, pEdu.pVessel.vMean);
Mxy_flat = getMxy_ss(Mz_flat, pEdu.pMri, pEdu.pMri.relax.blood);


pVesselLami          = pEdu.pVessel;
pVesselLami.profile  = 'parabolic1';
pVesselLami.S.lumen    = Mxy_flat;
pVesselLami.S.surround = Mxy_flat/nnz(pVesselLami.mask.lumen)*nnz(pVesselLami.mask.surround)/100;

pVesselPlug          = pEdu.pVessel;
pVesselPlug.profile  = 'plug';
pVesselPlug.S.lumen  = Mxy_flat;
pVesselPlug.S.surround = Mxy_flat/nnz(pVesselLami.mask.lumen)*nnz(pVesselLami.mask.surround)/100;


% Run simulations (light=false to retain magMap/vMap)
resPlug = runSim(pVesselPlug, pEdu.pSim, pEdu.pMri, [], false);
resLami = runSim(pVesselLami, pEdu.pSim, pEdu.pMri, [], false);


% Layout: 4 rows × 3 cols, taller figure so tiles approach square
%   col 1      : images (vel row 1|3, mag row 2|4)
%   cols 2-3   : complex-domain [2×2] (plug rows 1-2, lami rows 3-4)
fSim = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
hTSim = tiledlayout(fSim,4,3,'TileSpacing','compact','Padding','compact'); ax = {};
coorPE = resPlug.pSim.spinGrid.coorPE;
coorFE = resPlug.pSim.spinGrid.coorFE;
cVes   = [pVesselPlug.ID/2 * cos(theta); pVesselPlug.ID/2 * sin(theta)];

% Plug vel — row 1, col 1
ax{end+1} = nexttile(1);
imagesc(ax{end}, coorPE, coorFE, resPlug.vMap); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end},'Location','westoutside'), 'velocity [cm/s]');
title(ax{end}, 'plug | vel');

% Plug mag — row 2, col 1
ax{end+1} = nexttile(4);
imagesc(ax{end}, coorPE, coorFE, resPlug.magMap); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end},'Location','westoutside'), 'MR magn. [a.u.]');
title(ax{end}, 'plug | mag');

% Plug complex domain — rows 1-2, cols 2-3
ax{end+1} = nexttile(2,[2 2]);
plotComplexDomain(ax{end}, resPlug.I, resPlug.pMri.venc.vencList, 'full', 'line');
hold(ax{end},'on');
plot(ax{end},real(mean(resPlug.Is)),imag(mean(resPlug.Is)),'.r')
set(ax{end},'Box','on','XTick',[],'YTick',[]);
title(ax{end}, 'plug flow');

% Lami vel — row 3, col 1
ax{end+1} = nexttile(7);
imagesc(ax{end}, coorPE, coorFE, resLami.vMap); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end},'Location','westoutside'), 'velocity [cm/s]');
title(ax{end}, 'laminar | vel');

% Lami mag — row 4, col 1
ax{end+1} = nexttile(10);
imagesc(ax{end}, coorPE, coorFE, resLami.magMap); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end},'Location','westoutside'), 'MR magn. [a.u.]');
title(ax{end}, 'laminar | mag');

% Lami complex domain — rows 3-4, cols 2-3
ax{end+1} = nexttile(8,[2 2]);
plotComplexDomain(ax{end}, resLami.I, resLami.pMri.venc.vencList, 'full', 'line');
hold(ax{end},'on');
plot(ax{end},real(mean(resLami.Is)),imag(mean(resLami.Is)),'.r')
set(ax{end},'Box','on','XTick',[],'YTick',[]);
title(ax{end}, 'laminar flow');

% ax{1,4}=vel, ax{2,5}=mag, ax{3,6}=complex
cLim = get([ax{[1 4]}],'CLim'); cLim = [-1 1].*max(abs([cLim{:}])); set([ax{[1 4]}],'CLim',cLim);
cLim = get([ax{[2 5]}],'CLim'); cLim = [ 0 1].*max(    [cLim{:}] ); set([ax{[2 5]}],'CLim',cLim);



set(findall(fSim,'Type','axes'),'FontSize',14);
set(findall(fSim,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec3fig,'simSummary.fig'),'file')
    saveas(        fSim, fullfile(sec3fig,'simSummary.fig'));
    exportgraphics(fSim, fullfile(sec3fig,'simSummary.png'));
    exportgraphics(fSim, fullfile(sec3fig,'simSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%
end







if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 - radial profiles and fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec5fig = fullfile(info.project.figures, '5-radial-profiles-and-fits');
if ~exist(sec5fig,'dir'); mkdir(sec5fig); end

cFlow          = mean(data(:,:,dataVenc==inf),3);
cFlowSpoiled   = mean(data(:,:,dataVenc==min(dataVenc)),3);
cNoFlow        = mean(dataNoFlow(:,:,dataVenc==inf),3);
cNoFlowSpoiled = mean(dataNoFlow(:,:,dataVenc==min(dataVenc)),3);
vFlow   = phase2vel(angle(mean(data(      :,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
vNoFlow = phase2vel(angle(mean(dataNoFlow(:,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
% rFlow   = rGrid;
% rNoFlow = rGrid;


% Compute some masks
m = abs(mean(data(:,:,dataVenc==inf),3));
[maskBloodOnly, maskWallOnly, maskNonBloodOnly, maskTissueOnly, maskWallLowMag] = makeVesselMasks(FEgrid, PEgrid, FEspacing, PEspacing, ID, OD, m);


% % Fit of v(r) velocity function of radius
% R = ID/2;
% [velFit, velFitFixR] = fitVelProfile(rGrid(maskBloodOnly), vFlow(maskBloodOnly), R);


% Fit blood [v(r) and m(v(r)) jointly with free centre offset (dx,dy)] and Fit wall and tissue [peicewise plateau-slope-plateau]
R = ID/2;
[velJointFit, magJoinFit, velJointFit1D] = fitMagVelProfile(rGrid(maskBloodOnly), pGrid(maskBloodOnly), vFlow(maskBloodOnly), abs(cFlow(maskBloodOnly)), abs(cNoFlow(maskBloodOnly)), R, 'joint', 4, true);
rWO   = R + OD/2-ID/2;
dr    = (FEspacing + PEspacing)/2;
pLow  = mean(abs(cFlow(maskWallOnly)));
pHigh = mean(abs(cFlow(maskTissueOnly)));
magWallTissueFit = fitMagBoundaryProfile(rGrid(maskNonBloodOnly), pGrid(maskNonBloodOnly), [velJointFit.FEoffset velJointFit.PEoffset], abs(cFlow(maskNonBloodOnly)), rWO, dr, pLow, pHigh);

% Recompute masks->Refit->Recompute masks
[maskBloodOnly, maskWallOnly, maskNonBloodOnly, maskTissueOnly] = makeVesselMasks(FEgrid - velJointFit.FEoffset, PEgrid - velJointFit.PEoffset, FEspacing, PEspacing, velJointFit.R*2, velJointFit.R*2+magWallTissueFit.dr*2);
R = velJointFit.R;
[velJointFit, magJoinFit, velJointFit1D] = fitMagVelProfile(rGrid(maskBloodOnly), pGrid(maskBloodOnly), vFlow(maskBloodOnly), abs(cFlow(maskBloodOnly)), abs(cNoFlow(maskBloodOnly)), R, 'joint', 4, true);
rWO   = magWallTissueFit.rWO;
dr    = (FEspacing + PEspacing)/2;
pLow  = magWallTissueFit.pLow;
pHigh = magWallTissueFit.pHigh;
magWallTissueFit = fitMagBoundaryProfile(rGrid(maskNonBloodOnly), pGrid(maskNonBloodOnly), [velJointFit.FEoffset velJointFit.PEoffset], abs(cFlow(maskNonBloodOnly)), rWO, dr, pLow, pHigh);
[maskBloodOnly, maskWallOnly, maskNonBloodOnly, maskTissueOnly] = makeVesselMasks(FEgrid - velJointFit.FEoffset, PEgrid - velJointFit.PEoffset, FEspacing, PEspacing, velJointFit.R*2, velJointFit.R*2+magWallTissueFit.dr*2);
rGridOffset = sqrt((FEgrid - velJointFit.FEoffset).^2 + (PEgrid - velJointFit.PEoffset).^2);

% Fit of m(r) magnitude function of radius for the spoiled low-venc data
R = velJointFit.R;
B = abs(mean(cNoFlowSpoiled(maskBloodOnly)));
[magSpoilFit, magSpoilFitFixR, magSpoilFitFixB, magSpoilFitFixBR] = fitMagProfile(rGrid(maskBloodOnly), abs(cFlowSpoiled(maskBloodOnly)), R, B, 4);





% figures
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 33 16]);
hT = tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact'); ax = {};
cx = velJointFit.R * cos(theta) + velJointFit.PEoffset;
cy = velJointFit.R * sin(theta) + velJointFit.FEoffset;

% ── radial profile of magnitudes ─────────
ax{end+1} = nexttile(1,[4 4]);
plot(ax{end}, double(rGridOffset), double(abs(cFlow)), '.w');
hold(ax{end},'on');
plot(ax{end}, double(rGridOffset), double(abs(cFlowSpoiled)), '.m');
r = linspace(0,velJointFit.R,100);
plot(ax{end}, r, magJoinFit(velJointFit1D(r)),'-w');
plot(ax{end}, r, magSpoilFitFixBR(r),'-m');
xline(velJointFit1D.R,':'); xline(magWallTissueFit.rWO,':');
axis square tight; xlabel('radial position [mm]'); ylabel('MR mag');

% ── no-flow venc==inf maps ─────────
im = mean(dataNoFlow(:,:,dataVenc==inf),3);
ax{end+1} = nexttile(5);
imagesc(ax{end}, PEpos, FEpos, abs(im)); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'MR magn. [a.u.]');
title(ax{end}, 'no-flow | venc=\infty | mag');

ax{end+1} = nexttile(6);
imagesc(ax{end}, PEpos, FEpos, angle(im), [-pi pi]); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'phase [rad]');
title(ax{end}, 'no-flow | venc=\infty | phase');

% ── flow venc==inf maps ─────────
im = mean(data(:,:,dataVenc==inf),3);
ax{end+1} = nexttile(11);
imagesc(ax{end}, PEpos, FEpos, abs(im)); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'MR magn. [a.u.]');
title(ax{end}, 'flow | venc=\infty | mag');

ax{end+1} = nexttile(12);
imagesc(ax{end}, PEpos, FEpos, angle(im), [-pi pi]); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'phase [rad]');
title(ax{end}, 'flow | venc=\infty | phase');

% ── flow venc==bestVenc maps ─────────
im = mean(data(:,:,dataVenc==bestVenc),3);
ax{end+1} = nexttile(17);
imagesc(ax{end}, PEpos, FEpos, abs(im)); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'MR magn. [a.u.]');
title(ax{end}, ['flow | venc=' num2str(bestVenc) ' cm/s | mag']);

ax{end+1} = nexttile(18);
imagesc(ax{end}, PEpos, FEpos, angle(im), [-pi pi]); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'phase [rad]');
title(ax{end}, ['flow | venc=' num2str(bestVenc) ' cm/s | phase']);

% ── flow venc==min maps ─────────
im = mean(data(:,:,dataVenc==2),3);
ax{end+1} = nexttile(23);
imagesc(ax{end}, PEpos, FEpos, abs(im)); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'MR magn. [a.u.]');
title(ax{end}, 'flow | venc=min | mag');

ax{end+1} = nexttile(24);
imagesc(ax{end}, PEpos, FEpos, angle(im), [-pi pi]); axis image;
ax{end}.Colormap = blueBlackRed; set(ax{end},'XTick',[],'YTick',[]);
ylabel(colorbar(ax{end}), 'phase [rad]');
title(ax{end}, 'flow | venc=min | phase');

% harmonize magnitude CLims across all mag tiles (ax{2,4,6,8})
cLim = get([ax{[2 4 6 8]}],'CLim'); cLim = [0 1].*max([cLim{:}]); set([ax{[2 4 6 8]}],'CLim',cLim);
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec5fig,'magPhaseMaps.fig'),'file')
    saveas(        f, fullfile(sec5fig,'magPhaseMaps.fig'));
    exportgraphics(f, fullfile(sec5fig,'magPhaseMaps.png'));
    exportgraphics(f, fullfile(sec5fig,'magPhaseMaps.svg'));
end










f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
hT = tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact'); ax = {};
v_plt = linspace(0, max(vFlow(maskBloodOnly)), 200);

% ── Magnitude images ────────────────────────
% phantom
ax{end+1} = nexttile(1);
imagesc(PEpos,FEpos,abs(cFlow),[0 max(abs(cFlow(maskBloodOnly)))]); axis image;
colormap(ax{end},'gray'); colorbar;
title(ax{end},'measured mag');
% matched sim
ax{end+1} = nexttile(2);
imagesc(PEpos,FEpos,reshape(magJoinFit(velJointFit(rGrid(:),pGrid(:))),size(rGrid))); axis image;
colormap(ax{end},'gray'); colorbar;
title(ax{end},'fitted mag');

cLim = get([ax{end-1:end}],'CLim'); cLim=[0 1].*max([cLim{:}]); set([ax{end-1:end}],'CLim',cLim);


% ── velocity images ───────────────────────────────
% phantom
ax{end+1} = nexttile(4);
imagesc(PEpos,FEpos,vFlow,[-1 1].*max(vFlow(maskBloodOnly))); axis image;
ax{end}.Colormap = blueBlackRed; colorbar;
title(ax{end},'measured vel');
% matched sim
ax{end+1} = nexttile(5);
imagesc(PEpos,FEpos,reshape(velJointFit(rGrid(:),pGrid(:)),size(rGrid))); axis image;
ax{end}.Colormap = blueBlackRed; colorbar;
title(ax{end},'fitted vel');

cLim = get([ax{end-1:end}],'CLim'); cLim=[-1 1].*max(abs([cLim{:}])); set([ax{end-1:end}],'CLim',cLim);

% ── magnitude function of velocity ────────────────
% phantom data scatter + fit line
ax{end+1} = nexttile(7);
plot(ax{end}, double(vFlow(maskBloodOnly)), double(abs(cFlow(maskBloodOnly))), '.w');
hold(ax{end},'on');
plot(ax{end}, v_plt, magJoinFit(v_plt), 'y-', 'LineWidth', 1.5);
xlabel(ax{end},'velocity [cm/s]'); ylabel(ax{end},'magnitude [a.u.]');
set(ax{end},'Color','k','GridColor',[0.5 0.5 0.5]); grid(ax{end},'on');
title(ax{end},'measured vel vs mag');
axis square

% sim data scatter
ax{end+1} = nexttile(8);
v = velJointFit(rGrid(:),pGrid(:));
m = magJoinFit(v);
plot(ax{end}, double(v(v>0)), double(m(v>0)), '.w');
xlabel(ax{end},'velocity [cm/s]'); ylabel(ax{end},'magnitude [a.u.]');
set(ax{end},'Color','k','GridColor',[0.5 0.5 0.5]); grid(ax{end},'on');
title(ax{end},'sim vel vs mag');
axis square

xLim=get([ax{end-1:end}],'XLim'); xLim=[0 1].*max([xLim{:}]); set([ax{end-1:end}],'XLim',xLim);
yLim=get([ax{end-1:end}],'YLim'); yLim=[0 1].*max([yLim{:}]); set([ax{end-1:end}],'YLim',yLim);




% ── radial profile of magnitudes ─────────
ax{end+1} = nexttile(3);
% phantom data
plot(ax{end}, double(rGridOffset), double(abs(cFlow)), '.w');
hold(ax{end},'on');
% blood fit
r = linspace(0,velJointFit.R,100);
plot(ax{end}, r, magJoinFit(velJointFit1D(r)),'-');
xline(velJointFit1D.R)
% wall-tissue fit
r = linspace(velJointFit.R,max(rGridOffset(:)),100);
plot(ax{end}, r, magWallTissueFit(r),'-');
xline(magWallTissueFit.rWO);
axis square
xlabel('radial position [mm]')
ylabel('MR mag')




% ── radial profile of velocities ─────────
ax{end+1} = nexttile(6);
% phantom data
plot(ax{end}, double(rGridOffset(maskBloodOnly|maskTissueOnly)), double(vFlow(maskBloodOnly|maskTissueOnly)), '.w');
hold on
% blood fit
r = linspace(0,max(rGridOffset(:)),100);
plot(ax{end}, r, velJointFit1D(r),'-');
xline(velJointFit1D.R)
xline(magWallTissueFit.rWO)
axis square
xlabel('radial position [mm]')
ylabel('velocity [cm/s]')



set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec5fig,'radialProfilesFits.fig'),'file')
    saveas(        f, fullfile(sec5fig,'radialProfilesFits.fig'));
    exportgraphics(f, fullfile(sec5fig,'radialProfilesFits.png'));
    exportgraphics(f, fullfile(sec5fig,'radialProfilesFits.svg'));
end
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec5fig,'radialProfilesFits.mat'),'file')
    save(fullfile(sec5fig,'radialProfilesFits.mat'), 'velJointFit', 'magJoinFit', 'velJointFit1D', 'magWallTissueFit');
end

% % Comparison figure: no-offset vs free-offset joint fit (requires velJointFit_noOff)
% thetaC       = linspace(0,2*pi,360);
% vel_noOff    = reshape(velJointFit_noOff(rGrid(:)),          size(rGrid));
% vel_withOff  = reshape(velJointFit(      rGrid(:),pGrid(:)), size(rGrid));
% mag_noOff    = reshape(magJoinFit_noOff( vel_noOff(:)),      size(rGrid));
% mag_withOff  = reshape(magJoinFit(       vel_withOff(:)),    size(rGrid));
% ... (see git history for full comparison figure code)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % section 5




if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 - phantom and matched simulation summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec4fig = fullfile(info.project.figures, '4-phantom-and-matched-simulation-summary');
if ~exist(sec4fig,'dir'); mkdir(sec4fig); end
sec5fig = fullfile(info.project.figures, '5-radial-profiles-and-fits');

% --- Load joint fits ---
load(fullfile(sec5fig,'radialProfilesFits.mat'), 'velJointFit', 'magJoinFit', 'velJointFit1D', 'magWallTissueFit');
m = abs(mean(data(:,:,dataVenc==inf),3));
[maskBloodOnly, ~, ~, maskTissueOnly, maskWallLowMag] = makeVesselMasks(FEgrid - velJointFit.FEoffset, PEgrid - velJointFit.PEoffset, FEspacing, PEspacing, velJointFit.R*2, velJointFit.R*2+magWallTissueFit.dr*2, m);

% --- Simulation: setup matched to joint fit ---
p       = runSim;
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

pSim.voxGrid.fovFE = size(data,1) * FEspacing;
pSim.voxGrid.fovPE = size(data,2) * PEspacing;
pSim.voxGrid.matFE = size(data,1);
pSim.voxGrid.matPE = size(data,2);
pSim.nSpin         = (2^10)^2;
pSim.gridMode      = 'pseudoVoxel';

pVessel.ID      = velJointFit.R * 2;
pVessel.WT      = magWallTissueFit.rWO - velJointFit.R;
pVessel.vMean   = velJointFit.Vmax / 2;
pVessel.posFE   = velJointFit.FEoffset;
pVessel.posPE   = velJointFit.PEoffset;
pVessel.profile = 'parabolic1';

pMri.fieldStrength  = 3;
pMri.species        = 'phantom';
pMri.sliceThickness = 2.2;
pMri.TR             = 75.90/(5+1)/1000;
pMri.TE             = 9.8/1000;
pMri.FA             = 50;
pMri.venc.method    = 'FVEmono';
pMri.venc.FVEbw = 100;

p       = runSim(pVessel, pSim, pMri);
pSim    = p.pSim;
pVessel = p.pVessel;
pMri    = p.pMri;

% reconstruct cartesian and polar representation grids
[pSim.spinGrid.feGrid, pSim.spinGrid.peGrid] = ndgrid(pSim.spinGrid.coorFE, pSim.spinGrid.coorPE);  % dim1=FE, dim2=PE
pSim.spinGrid.rGrid = sqrt(pSim.spinGrid.peGrid.^2+pSim.spinGrid.feGrid.^2);
pSim.spinGrid.pGrid = -atan2(pSim.spinGrid.feGrid, pSim.spinGrid.peGrid);  % theta=0 → +PE ("right" in imagesc display)
[pSim.voxGrid.feGrid, pSim.voxGrid.peGrid] = ndgrid(pSim.voxGrid.coorFE, pSim.voxGrid.coorPE);  % dim1=FE, dim2=PE
pSim.voxGrid.rGrid = sqrt(pSim.voxGrid.peGrid.^2+pSim.voxGrid.feGrid.^2);
pSim.voxGrid.pGrid = -atan2(pSim.voxGrid.feGrid, pSim.voxGrid.peGrid);  % theta=0 → +PE ("right" in imagesc display)

% Lumen: per-spin signal from joint m(v(r,p)) fit
% velJointFit(r,p): r=distance from image centre, p=-atan2(FEgrid,PEgrid); offsets applied internally
pVessel.S.lumen    = max(0, magJoinFit(velJointFit(pSim.spinGrid.rGrid(pVessel.mask.lumen), pSim.spinGrid.pGrid(pVessel.mask.lumen))));

% Wall and surround: per-spin signal from boundary profile fit (function of offset-corrected radius)
pVessel.S.wall     = 0;
pVessel.S.surround = magWallTissueFit.pHigh;

% Run matched simulation
resMatchedSim = runSim(pVessel, pSim, pMri, [], false);



%
% Figure 1: combined summary — 4×4 row-major
%   col 1: phantom maps  [2×1 each], colorbar west
%   col 2: sim maps      [2×1 each], no colorbar
%   cols 3-4: complex domain [4×2]
%   tiles square at 35×18 cm → [2×1] ≈ 8.75×9 cm, [4×2] ≈ 17.5×18 cm
%
fComb  = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 32 16]);
hTComb = tiledlayout(fComb, 4, 4, 'TileSpacing','compact','Padding','compact');
axComb = {};

% Phantom mag — rows 1-2, col 1
axComb{end+1} = nexttile(hTComb, 1, [2 1]);
im = abs(mean(data(:,:,dataVenc==inf),3));
imagesc(axComb{end}, PEpos, FEpos, im, [0 max(im(:))]); axis image;
ylabel(colorbar(axComb{end},'Location','westoutside'), 'MR magn. [a.u.]');
axComb{end}.Colormap = gray; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, 'phantom ROI');

% Phantom vel — rows 3-4, col 1
axComb{end+1} = nexttile(hTComb, 9, [2 1]);
im = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
imagesc(axComb{end}, PEpos, FEpos, im, [-1 1].* max(abs(im(:)))); axis image;
ylabel(colorbar(axComb{end},'Location','westoutside'), 'velocity [cm/s]');
axComb{end}.Colormap = blueBlackRed; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, ['venc=' num2str(bestVenc) ' cm/s']);

% Sim mag — rows 1-2, col 2 (no colorbar)
axComb{end+1} = nexttile(hTComb, 2, [2 1]);
im = resMatchedSim.magMap;
imagesc(axComb{end}, resMatchedSim.pSim.spinGrid.coorPE, resMatchedSim.pSim.spinGrid.coorFE, im, [0 max(im(:))]); axis image;
axComb{end}.Colormap = gray; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, 'simulation ROI');

% Sim vel — rows 3-4, col 2 (no colorbar)
axComb{end+1} = nexttile(hTComb, 10, [2 1]);
im = resMatchedSim.vMap;
imagesc(axComb{end}, resMatchedSim.pSim.spinGrid.coorPE, resMatchedSim.pSim.spinGrid.coorFE, im, [-1 1].*max(abs(im(:)))); axis image;
axComb{end}.Colormap = blueBlackRed; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, sprintf('joint fit: vMean=%.1f cm/s, R=%.1f mm', velJointFit.Vmax/2, velJointFit.R));

% Complex domain — rows 1-4, cols 3-4
axComb{end+1} = nexttile(hTComb, 3, [4 2]);
plotComplexDomain(axComb{end}, sum(data, [1 2]), dataVenc, 'full', 'markers');
hold(axComb{end}, 'on');
plot(axComb{end}, real(squeeze(resMatchedSim.I)), imag(squeeze(resMatchedSim.I)), 'c-', 'LineWidth', 1.5);
hold(axComb{end}, 'off');
legend(axComb{end}, {'phantom','simulation'}, 'Location','best', 'TextColor','w', 'Color','k');
title(axComb{end}, 'complex-domain ROI signal');

set(findall(fComb,'Type','axes'),'FontSize',14);
set(findall(fComb,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec4fig,'phantomSimCombinedSummary.fig'),'file')
    saveas(        fComb, fullfile(sec4fig,'phantomSimCombinedSummary.fig'));
    exportgraphics(fComb, fullfile(sec4fig,'phantomSimCombinedSummary.png'));
    exportgraphics(fComb, fullfile(sec4fig,'phantomSimCombinedSummary.svg'));
end

%
% Figure 2: masks — 2×5 rowmajor, no colorbars so all tiles have identical physical mm size
%  Row 1 (phantom): mag | vel | blood mask | tissue mask | wall mask
%  Row 2 (simulation): mag | vel | lumen mask | wall mask | surround mask
%  XLim/YLim forced equal within each row after axis image.
%
fMask  = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 16]);
hTMask = tiledlayout(fMask, 2, 5, 'TileSpacing','compact','Padding','compact');
axMask = {};
theta3 = linspace(0, 2*pi, 360);
PEax3  = resMatchedSim.pSim.spinGrid.coorPE;
FEax3  = resMatchedSim.pSim.spinGrid.coorFE;

% Phantom row — no colorbars
axMask{end+1} = nexttile(hTMask);
im = abs(mean(data(:,:,dataVenc==inf),3));
imagesc(axMask{end}, PEpos, FEpos, im, [0 max(im(:))]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, 'phantom mag');

axMask{end+1} = nexttile(hTMask);
im = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)), vencToM1(bestVenc));
imagesc(axMask{end}, PEpos, FEpos, im, [-1 1].*max(abs(im(:)))); axis image;
axMask{end}.Colormap = blueBlackRed; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, ['phantom vel  venc=' num2str(bestVenc) ' cm/s']);

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskBloodOnly, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
title(axMask{end}, 'blood mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskTissueOnly, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
title(axMask{end}, 'tissue mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskWallLowMag, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
hold(axMask{end},'on');
title(axMask{end}, 'wall mask');

% Simulation row — no colorbars
axMask{end+1} = nexttile(hTMask);
im = resMatchedSim.magMap;
imagesc(axMask{end}, PEax3, FEax3, im, [0 max(im(:))]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, 'simulation mag');

axMask{end+1} = nexttile(hTMask);
im = resMatchedSim.vMap;
imagesc(axMask{end}, PEax3, FEax3, im, [-1 1].*max(abs(im(:)))); axis image;
axMask{end}.Colormap = blueBlackRed; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, sprintf('simulation vel  vMean=%.1f cm/s', velJointFit.Vmax/2));

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.lumen), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, 'lumen mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.wall), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
hold(axMask{end},'on');
title(axMask{end}, 'wall mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.surround), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, 'surround mask');

% Force identical XLim/YLim within each row so all tiles show exactly the same mm extent
% (no colorbars means equal tile area; matching limits guarantees equal mm scale)
phXLim = axMask{1}.XLim; phYLim = axMask{1}.YLim;
for k = 2:5; set(axMask{k}, 'XLim',phXLim, 'YLim',phYLim); end
simXLim = axMask{6}.XLim; simYLim = axMask{6}.YLim;
for k = 7:10; set(axMask{k}, 'XLim',simXLim, 'YLim',simYLim); end

set(findall(fMask,'Type','axes'),'FontSize',14);
set(findall(fMask,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec4fig,'phantomSimMasksSummary.fig'),'file')
    saveas(        fMask, fullfile(sec4fig,'phantomSimMasksSummary.fig'));
    exportgraphics(fMask, fullfile(sec4fig,'phantomSimMasksSummary.png'));
    exportgraphics(fMask, fullfile(sec4fig,'phantomSimMasksSummary.svg'));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % section 4




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in vivo data -- sub-01 and sub-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%
%% 6 - in vivo summary
%%%%%%%%%%%%%%%%%%%%%%
sec6fig = fullfile(info.project.figures, '6-in-vivo-summary');
if ~exist(sec6fig,'dir'); mkdir(sec6fig); end

% Choose a subject and vessel
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
% vencSel_iv  = vencList_iv(round(end/2));
vencSel_iv  = 10;
refAvg_iv   = squeeze(mean(img(:,:,:,:,:,:,imgInfo.vencList==inf,         :,:,:,:,:,:,:,:,:), [7 11]));
selAvg_iv   = squeeze(mean(img(:,:,:,:,:,:,imgInfo.vencList==vencSel_iv,  :,:,:,:,:,:,:,:,:), [7 11]));
velMap_iv   = phase2vel(angle(selAvg_iv ./ refAvg_iv), vencToM1(vencSel_iv));

% ROI overlay -- separate figure
hFroi = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
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
set(findall(hFroi,'Type','axes'),'FontSize',14);
set(findall(hFroi,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec6fig,[figName_roi '.png']),'file')
    drawnow;
    exportgraphics(hFroi, fullfile(sec6fig,[figName_roi '.png']));
    exportgraphics(hFroi, fullfile(sec6fig,[figName_roi '.svg']));
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
% Layout: 6 rows × 10 cols
%   left  6×4: full-slice [4×4 top], ROI mag [2×2 bot-left], vel [2×2 bot-right]
%   right 6×6: complex-domain plot
f_ivs = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 27 16]);
hT_ivs = tiledlayout(f_ivs,6,10,'TileSpacing','compact','Padding','compact');
ax_ivs = {};

% Full-slice reference magnitude with ROI box — top 4×4
ax_ivs{end+1} = nexttile(hT_ivs, 1, [4 4]);
imagesc(ax_ivs{end}, refMag_iv, [0, 0.5*max(refMag_iv(:))]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
hold(ax_ivs{end},'on');
plot(ax_ivs{end}, [roiX(1)-.5 roiX(2)+.5 roiX(2)+.5 roiX(1)-.5 roiX(1)-.5], ...
                  [roiY(1)-.5 roiY(1)-.5 roiY(2)+.5 roiY(2)+.5 roiY(1)-.5], 'c', 'LineWidth', 1);
title(ax_ivs{end}, [subName ' vessel 01']);

% ROI crop -- reference magnitude — bottom-left 2×2
ax_ivs{end+1} = nexttile(hT_ivs, 41, [2 2]);
refMagCrop_iv = refMag_iv(roiY(1):roiY(2), roiX(1):roiX(2));
imagesc(ax_ivs{end}, refMagCrop_iv, [0, max(refMagCrop_iv(:))]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
xlabel(colorbar(ax_ivs{end},'Location','southoutside'), 'MR magn. [a.u.]');

% ROI crop -- velocity map — bottom-right 2×2
ax_ivs{end+1} = nexttile(hT_ivs, 43, [2 2]);
velCrop_iv = velMap_iv(roiY(1):roiY(2), roiX(1):roiX(2));
vLim_ivs = max(abs(velCrop_iv(:)));
imagesc(ax_ivs{end}, velCrop_iv, [-vLim_ivs vLim_ivs]);
set(ax_ivs{end},'XTick',[],'YTick',[],'Colormap',blueBlackRed,'DataAspectRatio',[imgInfo.res 1]);
xlabel(colorbar(ax_ivs{end},'Location','southoutside'), 'velocity [cm/s]');

% Complex-domain signal — right 6×6
ax_ivs{end+1} = nexttile(hT_ivs, 5, [6 6]);
plotComplexDomain(ax_ivs{end}, trjIV(:), [], 'full', 'markers');
set(ax_ivs{end}, 'XTick', [], 'YTick', [], 'Box', 'on');

% Save
figName_ivs = [subName '-vessel01-summary'];
set(findall(f_ivs,'Type','axes'),'FontSize',14);
set(findall(f_ivs,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec6fig,[figName_ivs '.fig']),'file')
    saveas(        f_ivs, fullfile(sec6fig,[figName_ivs '.fig']));
    exportgraphics(f_ivs, fullfile(sec6fig,[figName_ivs '.png']));
    exportgraphics(f_ivs, fullfile(sec6fig,[figName_ivs '.svg']));
end



%% %%%%%%%%%%%%%%%%%%%
end % section 6




if 1
%%%%%%%%%%%%%%%%%
%% 7 - conclusion
%%%%%%%%%%%%%%%%%%
% No MATLAB output — poster section only.
%% %%%%%%%%%%%%%%
end




if 1
saveThis = 1;
%%%%%%%%%%
%% 8 - FVE
%%%%%%%%%%
sec8fig = fullfile(info.project.figures, '8-extra');
if ~exist(sec8fig,'dir'); mkdir(sec8fig); end
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

fVelSpec = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
hFlat = plot(resSatin.pMri.venc.FVEvel,abs(fftshift(fft(squeeze(resSatin.I)))),'w');
hold on
hVdep = plot(res.pMri.venc.FVEvel     ,abs(fftshift(fft(squeeze(res.I))))     ,'g');
axis tight; grid on; xlabel('velocity (cm/s)'); ylabel('velocity spectrum/histogram');
yLim = ylim; yLim(1) = 0; ylim(yLim);
[N,edges] = histcounts(res.vMap(getVoxIdx(res.pSim.voxGrid,res.pSim.spinGrid)==0),20);
hVhist = histogram('BinEdges',edges,'BinCounts',N/max(N)*yLim(2),'FaceColor',0.5.*[1 1 1],'EdgeColor','none');
legend([hFlat,hVdep,hVhist],['velocity spectrum from' newline 'flat magnitude profile'],['velocity spectrum from' newline 'velocity-dependent magnitude profile'],'normalized velocity histogram','Location','northwest','box','off');
uistack(hVhist,'bottom');

fVelSpecInflow = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
% inflowVel = linspace(0,p.pVessel.vMean*3,2^10);
inflowVel = linspace(0,6,2^10);
[inflowMz,~,~,~,~,inflowVel] = getMz_ss(p.pMri,p.pMri.relax.blood,inflowVel);
% inflowMxy = getMxy_ss(inflowMz,p.pMri,p.pMri.relax.blood);
hStairs = stairs(inflowVel,inflowMz,'g');
axis tight square; grid on; xlabel('spin velocity (cm/s)'); ylabel('M_z');
ylim([0 1])

set(findall(fVelSpec,      'Type','axes'),'FontSize',14);
set(findall(fVelSpec,      'Type','text'),'FontSize',8);
set(findall(fVelSpecInflow,'Type','axes'),'FontSize',14);
set(findall(fVelSpecInflow,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec8fig,'FVEvelSpec.fig'),'file') || ~exist(fullfile(sec8fig,'FVEvelSpec_inflow.fig'),'file')
    saveas(        fVelSpec      , fullfile(sec8fig,'FVEvelSpec.fig'       ));
    exportgraphics(fVelSpec      , fullfile(sec8fig,'FVEvelSpec.png'       ));
    exportgraphics(fVelSpec      , fullfile(sec8fig,'FVEvelSpec.svg'       ));
    saveas(        fVelSpecInflow, fullfile(sec8fig,'FVEvelSpec_inflow.fig'));
    exportgraphics(fVelSpecInflow, fullfile(sec8fig,'FVEvelSpec_inflow.png'));
    exportgraphics(fVelSpecInflow, fullfile(sec8fig,'FVEvelSpec_inflow.svg'));
end
% FVE spectra reflects spin velocity distribution, but weighted by velocity-dependent spin magnitude
fVelSpec;
% Here the weighting effect was maximized using a 90 flip angle for a linear magnitude function of velocity
fVelSpecInflow;
%% %%%%%%%
end % section 8




if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9 - Phantom details -- all maps and masks with and without flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec9fig = fullfile(info.project.figures, '9-phantom-details');
if ~exist(sec9fig,'dir'); mkdir(sec9fig); end
m = abs(mean(data(:,:,dataVenc==inf),3));
[maskBloodOnly, ~, ~, maskTissueOnly, maskWallLowMag] = makeVesselMasks(FEgrid, PEgrid, FEspacing, PEspacing, ID, OD, m);
% Mag flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'magFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'magFlowOn.fig'       ));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOn.png'       ));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOn.svg'       ));
end



% Mag flow off
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'magFlowOff.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'magFlowOff.fig'));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOff.png'));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOff.svg'));
end


% Phase flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
        cb.Ticks = -pi:pi/2:pi;
        cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    end
    ax{end}.Colormap = blueBlackRed;
end

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'phaseFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'phaseFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOn.svg'));
end



% Phase flow off
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
        cb.Ticks = -pi:pi/2:pi;
        cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    end
    ax{end}.Colormap = blueBlackRed;
end

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'phaseFlowOff.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'phaseFlowOff.fig'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOff.png'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOff.svg'));
end




% CDvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
    ax{end}.Colormap = blueBlackRed;
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
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'CDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'CDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'CDvelFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'CDvelFlowOn.svg'));
end



% PDvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
    ax{end}.Colormap = blueBlackRed;
    % cb.Ticks = -pi:pi/2:pi;
    % cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
end
set([ax{:}],'CLim',[-1 1].*9);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'PDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'PDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'PDvelFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'PDvelFlowOn.svg'));
end




% Pvel flow on
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
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
    ax{end}.Colormap = blueBlackRed;
    % cb.Ticks = -pi:pi/2:pi;
    % cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
end
set([ax{:}],'CLim',[-1 1].*9);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly)
ax{end}.Colormap = gray;
hold(ax{end},'on');
title(ax{end},'blood-only mask');

%save
set(findall(f,'Type','axes'),'FontSize',14);
set(findall(f,'Type','text'),'FontSize',8);
if saveThis || ~exist(fullfile(sec9fig,'PvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'PvelFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'PvelFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'PvelFlowOn.svg'));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % section 9



if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10 - in vivo details -- all vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec10fig = fullfile(info.project.figures, '10-in-vivo-details');
if ~exist(sec10fig,'dir'); mkdir(sec10fig); end
for s = 1:2
    img      = inVivoSubData{s}.img;
    imgInfo  = inVivoSubData{s}.imgInfo;
    refImgAv = inVivoSubData{s}.refImgAv;

    % ROI overlay figure
    hF = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
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
    set(findall(hF,'Type','axes'),'FontSize',14);
    set(findall(hF,'Type','text'),'FontSize',8);
    if saveThis || ~exist(fullfile(sec10fig,[inVivoSubNames{s} '-roiOverlay.png']),'file')
        drawnow;
        exportgraphics(hF, fullfile(sec10fig,[inVivoSubNames{s} '-roiOverlay.png']));
        exportgraphics(hF, fullfile(sec10fig,[inVivoSubNames{s} '-roiOverlay.svg']));
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

        hFv = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
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
        thetaUnit = linspace(0,2*pi,100);
        x = abs(mean(trj(trjVenc==inf,:),[1 2]))*cos(thetaUnit);
        y = abs(mean(trj(trjVenc==inf,:),[1 2]))*sin(thetaUnit);
        uistack(plot(x,y,'w'),'bottom');
        title(legend(hPcont,trjVencLabel),'V_{enc} (cm/s)');
        vesselFigName = [inVivoSubNames{s} '_vessel-' num2str(roiIdx,'%02d')];
        subFigDir = fullfile(sec10fig, inVivoSubNames{s});
        if ~exist(subFigDir,'dir'); mkdir(subFigDir); end
        set(findall(hFv,'Type','axes'),'FontSize',14);
        set(findall(hFv,'Type','text'),'FontSize',8);
        if saveThis || ~exist(fullfile(subFigDir,[vesselFigName '.png']),'file')
            drawnow;
            exportgraphics(hFv, fullfile(subFigDir,[vesselFigName '.png']));
            exportgraphics(hFv, fullfile(subFigDir,[vesselFigName '.svg']));
        end
        close(hFv); clear hPcont d prob
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % section 10
