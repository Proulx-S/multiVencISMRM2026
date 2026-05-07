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

% Compute some masks
dFE = abs(FEgrid);
dPE = abs(PEgrid);
% farthest corner of each pixel from the center
d_far  = sqrt((dFE + FEspacing/2).^2 + (dPE + PEspacing/2).^2);
% nearest point of each pixel to the center
d_near = sqrt(max(0, dFE - FEspacing/2).^2 + max(0, dPE - PEspacing/2).^2);

maskBloodOnly    = d_far  < ID/2;              % pixel entirely inside inner circle
maskWallOnly     = d_near > ID/2 & d_far < OD/2; % pixel entirely within wall annulus
maskNonBloodOnly = d_near > ID/2;
maskTissueOnly   = d_near > OD/2;              % pixel entirely outside outer circle
maskWallLowMag   = M<min(M(maskTissueOnly));          % low magnitude pixels

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


% Flat magnitude profile: S.lumen = constant Mxy at vMean (no velocity-dependent inflow saturation)
Mz_flat  = getMz_ss(pEdu.pMri, pEdu.pMri.relax.blood, pEdu.pVessel.vMean);
Mxy_flat = getMxy_ss(Mz_flat, pEdu.pMri, pEdu.pMri.relax.blood);

pVesselPara          = pEdu.pVessel;
pVesselPara.S.lumen  = Mxy_flat;

pVesselPlug          = pEdu.pVessel;
pVesselPlug.PD       = pVesselPlug.ID;
pVesselPlug.S.lumen  = Mxy_flat;

% PCmono with venc list matching in vivo data range
pMri = pEdu.pMri;
pMri.venc.method  = 'PCmono';

% Run simulations (light=false to retain magMap/vMap)
resPara = runSim(pVesselPara, pEdu.pSim, pMri, [], false);
resPlug = runSim(pVesselPlug, pEdu.pSim, pMri, [], false);

% Build multi-venc signal arrays for plotMultiVenc
% res.I dims: [1 1 1 1 nVenc 2] = [FE PE SL t M1 M1ref]
% col 1 = velocity-encoded, col 2 = M1=0 reference (real after phase subtraction)
vencList = pMri.venc.vencList;
resList  = {resPlug, resPara};
flowNames = {'plug flow','laminar flow'};
Iplot    = cell(1,2);
vPlot    = cell(1,2);
for flowIdx = 1:2
    r     = resList{flowIdx};
    Iref  = r.I(1,1,1,1,1,2);
    Ienc  = squeeze(r.I(1,1,1,1,:,1));
    Iplot{flowIdx} = [Iref; Ienc(:)];
    vPlot{flowIdx} = [inf;  vencList];
end

% Figure: 2 rows (plug | laminar) x 3 cols (magMap | vMap | complex-plane spiral)
fSim = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 22]);
hTSim = tiledlayout(fSim,2,3,'TileSpacing','compact','Padding','compact'); axSim = {};
for rowIdx = 1:2
    r_ = resList{rowIdx};
    axSim{end+1} = nexttile(hTSim);
    imagesc(r_.magMap); axis image;
    axSim{end}.Colormap = gray;
    axSim{end}.CLim = [0 max(r_.magMap(:))*1.1];
    set(axSim{end},'XTick',[],'YTick',[]);
    title(axSim{end},{flowNames{rowIdx},'magnitude map'});
    axSim{end+1} = nexttile(hTSim);
    vLim = max(abs(r_.vMap(:)))*1.1;
    imagesc(r_.vMap,[-vLim vLim]); axis image;
    axSim{end}.Colormap = redblue;
    set(axSim{end},'XTick',[],'YTick',[]);
    ylabel(colorbar,'velocity (cm/s)','FontSize',8);
    title(axSim{end},{flowNames{rowIdx},'velocity map'});
    axSim{end+1} = nexttile(hTSim);
    I_s    = Iplot{rowIdx};
    v_s    = vPlot{rowIdx};
    Mnorm = abs(mean(I_s(v_s==inf)));
    plotComplexDomain(axSim{end}, I_s/Mnorm, [], 'tight', 'line');
    title(axSim{end},{flowNames{rowIdx},'complex-domain signal'});
end
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
vFlow   = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
vNoFlow = zeros(size(vFlow));
rFlow   = rGrid;
rNoFlow = rGrid;

% % Fit of v(r) velocity function of radius
% R = ID/2;
% [velFit, velFitFixR] = fitVelProfile(rGrid(maskBloodOnly), vFlow(maskBloodOnly), R);

% % Fit of m(r) magnitude function of radius
% R = ID/2;
% B = abs(mean(cNoFlow(maskBloodOnly)));
% [magFit, magFitFixR, magFitFixB, magFitFixBR] = fitMagProfile(rGrid(maskBloodOnly), abs(cFlow(maskBloodOnly)), R, B);

% % Fit of m(r) magnitude function of radius for the spoiled low-venc data
% R = ID/2;
% B = abs(mean(cNoFlowSpoiled(maskBloodOnly)));
% [magSpoilFit, magSpoilFitFixR, magSpoilFitFixB, magSpoilFitFixBR] = fitMagProfile(rGrid(maskBloodOnly), abs(cFlowSpoiled(maskBloodOnly)), R, B);

% % Fit of v(r) and m(v(r)) jointly — no offset (comparison baseline)
% R = ID/2;
% [velJointFit_noOff, magJoinFit_noOff] = fitMagVelProfile(rGrid(maskBloodOnly), vFlow(maskBloodOnly), abs(cFlow(maskBloodOnly)), abs(cNoFlow(maskBloodOnly)), R, 'joint', 4);

% Fit of v(r) and m(v(r)) jointly — free centre offset (dx,dy)
R = ID/2;
[velJointFit, magJoinFit, velJointFit1D] = fitMagVelProfile(rGrid(maskBloodOnly), pGrid(maskBloodOnly), vFlow(maskBloodOnly), abs(cFlow(maskBloodOnly)), abs(cNoFlow(maskBloodOnly)), R, 'joint', 4, true);

% Piecewise fit of amplitude radial profile (non-blood voxels) to find rWO.
% Model: low plateau | linear slope centred at rWO over width dr | high plateau
%   th = [pLow, pHigh, rWO, dr]
rWO   = R + OD/2-ID/2;
dr    = (FEspacing + PEspacing)/2;
pLow  = mean(abs(cFlow(maskWallOnly)));
pHigh = mean(abs(cFlow(maskTissueOnly)));
magWallTissueFit = fitMagBoundaryProfile(rGrid(maskNonBloodOnly), pGrid(maskNonBloodOnly), [velJointFit.FEoffset velJointFit.PEoffset], abs(cFlow(maskNonBloodOnly)), rWO, dr, pLow, pHigh);



% Offset-corrected radius for all pixels — used for radial profile plots
r_off_all = sqrt((FEgrid - velJointFit.FEoffset).^2 + (PEgrid - velJointFit.PEoffset).^2);

f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 30]);
hT = tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact'); ax = {};
theta  = linspace(0, 2*pi, 360);
cx     = velJointFit.R * cos(theta) + velJointFit.PEoffset;  % fitted lumen circle (imagesc: x=PE,y=FE)
cy     = velJointFit.R * sin(theta) + velJointFit.FEoffset;
v_plt  = linspace(0, max(vFlow(maskBloodOnly)), 200);
vBlood = velJointFit(rGrid(maskBloodOnly), pGrid(maskBloodOnly));

r_blood     = linspace(0, velJointFit.R, 200);
m_blood_fit = magJoinFit(velJointFit1D(r_blood));
r_all_plt   = linspace(0, max(r_off_all(:)), 400);

% ── Magnitude images ────────────────────────
% phantom
ax{end+1} = nexttile(1);
imagesc(PEpos,FEpos,abs(cFlow),[0 max(abs(cFlow(maskBloodOnly)))]); axis image;
colormap(ax{end},'gray'); colorbar; hold on; plot(cx,cy,'w--','LineWidth',1);
title(ax{end},'measured mag');
% matched sim
ax{end+1} = nexttile(2);
imagesc(PEpos,FEpos,reshape(magJoinFit(velJointFit(rGrid(:),pGrid(:))),size(rGrid))); axis image;
colormap(ax{end},'gray'); colorbar; hold on; plot(cx,cy,'w--','LineWidth',1);
title(ax{end},'fitted mag');

cLim = get([ax{end-1:end}],'CLim'); cLim=[0 1].*max([cLim{:}]); set([ax{end-1:end}],'CLim',cLim);


% ── velocity images ───────────────────────────────
% phantom
ax{end+1} = nexttile(4);
imagesc(PEpos,FEpos,vFlow,[-1 1].*max(vFlow(maskBloodOnly))); axis image;
colormap(ax{end},'jet'); colorbar; hold on; plot(cx,cy,'w--','LineWidth',1);
title(ax{end},'measured vel');
% matched sim
ax{end+1} = nexttile(5);
imagesc(PEpos,FEpos,reshape(velJointFit(rGrid(:),pGrid(:)),size(rGrid))); axis image;
colormap(ax{end},'jet'); colorbar; hold on; plot(cx,cy,'w--','LineWidth',1);
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
plot(ax{end}, double(r_off_all), double(abs(cFlow)), '.w');
hold(ax{end},'on');
% blood fit
r = linspace(0,velJointFit.R,100);
plot(ax{end}, r, magJoinFit(velJointFit1D(r)),'-');
xline(velJointFit1D.R)
% wall-tissue fit
r = linspace(velJointFit.R,max(r_off_all(:)),100);
plot(ax{end}, r, magWallTissueFit(r),'-');
xline(magWallTissueFit.rWO);
axis square
xlabel('radial position [mm]')
ylabel('MR mag')




% ── radial profile of velocities ─────────
ax{end+1} = nexttile(6);
% phantom data
plot(ax{end}, double(r_off_all(maskBloodOnly|maskTissueOnly)), double(vFlow(maskBloodOnly|maskTissueOnly)), '.w');
hold on
% blood fit
r = linspace(0,max(r_off_all(:)),100);
plot(ax{end}, r, velJointFit1D(r),'-');
xline(velJointFit1D.R)
xline(magWallTissueFit.rWO)
axis square
xlabel('radial position [mm]')
ylabel('velocity [cm/s]')



if saveThis || ~exist(fullfile(sec5fig,'radialProfilesFits.fig'),'file')
    saveas(        f, fullfile(sec5fig,'radialProfilesFits.fig'));
    exportgraphics(f, fullfile(sec5fig,'radialProfilesFits.png'));
    exportgraphics(f, fullfile(sec5fig,'radialProfilesFits.svg'));
end
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

% Offset-corrected radius for each spin (same convention as magWallTissueFit input)
r_off_spin = sqrt((pSim.spinGrid.feGrid - pVessel.posFE).^2 + (pSim.spinGrid.peGrid - pVessel.posPE).^2);

% Lumen: per-spin signal from joint m(v(r,p)) fit
% velJointFit(r,p): r=distance from image centre, p=-atan2(FEgrid,PEgrid); offsets applied internally
pVessel.S.lumen    = max(0, magJoinFit(velJointFit(pSim.spinGrid.rGrid(pVessel.mask.lumen), pSim.spinGrid.pGrid(pVessel.mask.lumen))));

% Wall and surround: per-spin signal from boundary profile fit (function of offset-corrected radius)
pVessel.S.wall     = 0;
pVessel.S.surround = magWallTissueFit.pHigh;

% Run matched simulation
resMatchedSim = runSim(pVessel, pSim, pMri, [], false);



%
% Figure 1: combined summary (no masks) — 2×5 columnmajor
%
fComb  = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 13]);
hTComb = tiledlayout(fComb, 2, 5, 'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor');
axComb = {};

% Phantom — column 1: mag, velocity
axComb{end+1} = nexttile(hTComb);
im = abs(mean(data(:,:,dataVenc==inf),3));
imagesc(axComb{end}, PEpos, FEpos, im, [0 max(im(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axComb{end}.Colormap = gray; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, 'phantom ROI');

axComb{end+1} = nexttile(hTComb);
im = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)),vencToM1(bestVenc));
imagesc(axComb{end}, PEpos, FEpos, im, [-1 1].* max(abs(im(:)))); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axComb{end}.Colormap = redblue; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, ['venc=' num2str(bestVenc) ' cm/s']);

% Simulation — column 2: mag, velocity
axComb{end+1} = nexttile(hTComb);
im = resMatchedSim.magMap;
imagesc(axComb{end}, resMatchedSim.pSim.spinGrid.coorPE, resMatchedSim.pSim.spinGrid.coorFE, im, [0 max(im(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
axComb{end}.Colormap = gray; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, 'simulation ROI');

axComb{end+1} = nexttile(hTComb);
im = resMatchedSim.vMap;
imagesc(axComb{end}, resMatchedSim.pSim.spinGrid.coorPE, resMatchedSim.pSim.spinGrid.coorFE, im, [-1 1].*max(abs(im(:)))); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
axComb{end}.Colormap = redblue; set(axComb{end},'XTick',[],'YTick',[]);
title(axComb{end}, sprintf('joint fit: vMean=%.1f cm/s, R=%.1f mm', velJointFit.Vmax/2, velJointFit.R));


% Combined complex-domain (cols 3-5): phantom dots + simulation line
axComb{end+1} = nexttile(hTComb, [2 3]);

plotComplexDomain(axComb{end}, sum(data, [1 2]), dataVenc, 'full', 'markers');
hold(axComb{end}, 'on');
plot(axComb{end}, real(squeeze(resMatchedSim.I)), imag(squeeze(resMatchedSim.I)), 'c-', 'LineWidth', 1.5);
hold(axComb{end}, 'off');
legend(axComb{end}, {'phantom','simulation'}, 'Location','best', 'TextColor','w', 'Color','k');
title(axComb{end}, 'complex-domain ROI signal');

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
fMask  = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 14]);
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
axMask{end}.Colormap = redblue; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, ['phantom vel  venc=' num2str(bestVenc) ' cm/s']);

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskBloodOnly, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
hold(axMask{end},'on'); plot(axMask{end}, ID/2*cos(theta3), ID/2*sin(theta3), 'm');
title(axMask{end}, 'blood mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskTissueOnly, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
hold(axMask{end},'on'); plot(axMask{end}, OD/2*cos(theta3), OD/2*sin(theta3), 'm');
title(axMask{end}, 'tissue mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEpos, FEpos, maskWallLowMag, [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[],'Color','none');
hold(axMask{end},'on');
plot(axMask{end}, ID/2*cos(theta3), ID/2*sin(theta3), 'm');
plot(axMask{end}, OD/2*cos(theta3), OD/2*sin(theta3), 'm');
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
axMask{end}.Colormap = redblue; set(axMask{end},'XTick',[],'YTick',[]);
title(axMask{end}, sprintf('simulation vel  vMean=%.1f cm/s', velJointFit.Vmax/2));

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.lumen), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
hold(axMask{end},'on'); plot(axMask{end}, ID/2*cos(theta3), ID/2*sin(theta3), 'm');
title(axMask{end}, 'lumen mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.wall), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
hold(axMask{end},'on');
plot(axMask{end}, ID/2*cos(theta3), ID/2*sin(theta3), 'm');
plot(axMask{end}, OD/2*cos(theta3), OD/2*sin(theta3), 'm');
title(axMask{end}, 'wall mask');

axMask{end+1} = nexttile(hTMask);
imagesc(axMask{end}, PEax3, FEax3, single(resMatchedSim.pVessel.mask.surround), [0 1]); axis image;
axMask{end}.Colormap = gray; set(axMask{end},'XTick',[],'YTick',[]);
hold(axMask{end},'on'); plot(axMask{end}, OD/2*cos(theta3), OD/2*sin(theta3), 'm');
title(axMask{end}, 'surround mask');

% Force identical XLim/YLim within each row so all tiles show exactly the same mm extent
% (no colorbars means equal tile area; matching limits guarantees equal mm scale)
phXLim = axMask{1}.XLim; phYLim = axMask{1}.YLim;
for k = 2:5; set(axMask{k}, 'XLim',phXLim, 'YLim',phYLim); end
simXLim = axMask{6}.XLim; simYLim = axMask{6}.YLim;
for k = 7:10; set(axMask{k}, 'XLim',simXLim, 'YLim',simYLim); end

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
plotComplexDomain(ax_ivs{end}, trjIV(:), [], 'tight', 'markers');

% Save
figName_ivs = [subName '-vessel01-summary'];
if saveThis || ~exist(fullfile(sec6fig,[figName_ivs '.fig']),'file')
    saveas(        f_ivs, fullfile(sec6fig,[figName_ivs '.fig']));
    exportgraphics(f_ivs, fullfile(sec6fig,[figName_ivs '.png']));
    exportgraphics(f_ivs, fullfile(sec6fig,[figName_ivs '.svg']));
end



%% %%%%%%%%%%%%%%%%%%%
end % section 6




if 1
%%%%%%%%%%%%%%%%%%
%% 7 - conclusion
%%%%%%%%%%%%%%%%%%
% No MATLAB output — poster section only.
%% %%%%%%%%%%%%
end




if 1
saveThis = 1;
%%%%%%%%%%
%% 8 - FVE
%%%%%%%%%%%%
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

fVelSpec = figure;
hFlat = plot(resSatin.pMri.venc.FVEvel,abs(fftshift(fft(squeeze(resSatin.I)))),'w');
hold on
hVdep = plot(res.pMri.venc.FVEvel     ,abs(fftshift(fft(squeeze(res.I))))     ,'g');
axis tight; grid on; xlabel('velocity (cm/s)'); ylabel('velocity spectrum/histogram');
yLim = ylim; yLim(1) = 0; ylim(yLim);
[N,edges] = histcounts(res.vMap(getVoxIdx(res.pSim.voxGrid,res.pSim.spinGrid)==0),20);
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
if saveThis || ~exist(fullfile(sec9fig,'magFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'magFlowOn.fig'       ));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOn.png'       ));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOn.svg'       ));
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
if saveThis || ~exist(fullfile(sec9fig,'magFlowOff.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'magFlowOff.fig'));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOff.png'));
    exportgraphics(f      , fullfile(sec9fig,'magFlowOff.svg'));
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
if saveThis || ~exist(fullfile(sec9fig,'phaseFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'phaseFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOn.svg'));
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
if saveThis || ~exist(fullfile(sec9fig,'phaseFlowOff.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'phaseFlowOff.fig'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOff.png'));
    exportgraphics(f      , fullfile(sec9fig,'phaseFlowOff.svg'));
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
if saveThis || ~exist(fullfile(sec9fig,'CDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'CDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'CDvelFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'CDvelFlowOn.svg'));
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
if saveThis || ~exist(fullfile(sec9fig,'PDvelFlowOn.fig'),'file')
    saveas(        f      , fullfile(sec9fig,'PDvelFlowOn.fig'));
    exportgraphics(f      , fullfile(sec9fig,'PDvelFlowOn.png'));
    exportgraphics(f      , fullfile(sec9fig,'PDvelFlowOn.svg'));
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
%%%%%%%%%%%%%%%
sec10fig = fullfile(info.project.figures, '10-in-vivo-details');
if ~exist(sec10fig,'dir'); mkdir(sec10fig); end
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
        subFigDir = fullfile(sec10fig, inVivoSubNames{s});
        if ~exist(subFigDir,'dir'); mkdir(subFigDir); end
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
