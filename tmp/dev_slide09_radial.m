% dev_slide09_radial.m -- dev script for radial profile fitting (slide 9)
% Run from: cd /scratch/bass/projects/multiVencISMRM2026 && matlab -batch "run('tmp/dev_slide09_radial.m')"

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');
addpath('/scratch/bass/projects/multiVenc'); % for fit_one_parameter_parabola

projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

saveThis = 1;

%% Load phantom data (same as doIt.m data loading section)
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
data = conj(data); dataNoFlow = conj(dataNoFlow);

ID =  6.35; % mm
OD = 11.11; % mm
bestVenc = 9; % cm/s

% Coordinates (recompute same as doIt.m)
M0 = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, size(M0,1)*FEspacing-FEspacing/2, size(M0,1));
PEpos = linspace(PEspacing/2, size(M0,2)*PEspacing-PEspacing/2, size(M0,2));
[FEgrid, PEgrid] = ndgrid(FEpos, PEpos);
total = sum(M0(:));
com(1) = sum(FEgrid(:).*M0(:))/total;
com(2) = sum(PEgrid(:).*M0(:))/total;
FEgrid = FEgrid-com(1); FEpos = FEpos-com(1);
PEgrid = PEgrid-com(2); PEpos = PEpos-com(2);
rGrid = sqrt(PEgrid.^2+FEgrid.^2);
dFE = abs(FEgrid); dPE = abs(PEgrid);
d_far  = sqrt((dFE+FEspacing/2).^2+(dPE+PEspacing/2).^2);
d_near = sqrt(max(0,dFE-FEspacing/2).^2+max(0,dPE-PEspacing/2).^2);
maskBloodOnly  = d_far  < ID/2;
maskTissueOnly = d_near > OD/2;
maskWallLowMag = single(M0<0.44e-7);

%% Compute phase difference and velocity at bestVenc
M  = abs(mean(data(:,:,dataVenc==inf),3));
PD    = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))), 3));
m1bestVenc = vencToM1(bestVenc); % convert venc (cm/s) to M1 (T*s^2/m)
velPD = phase2vel(PD, m1bestVenc);
CD = mean(data(:,:,dataVenc==bestVenc),3) - mean(data(:,:,dataVenc==inf),3);

%% Selection masks for fitting
idxBlood = maskBloodOnly; % pixels entirely inside vessel lumen

%% --- Parabolic fit to velocity vs. radius ---
rFit  = double(rGrid(idxBlood));  % [mm]
vFit  = double(velPD(idxBlood));  % [cm/s]
% Sign convention: after conj(), flow may be negative; work with the dominant sign
if mean(vFit) < 0; vFit = -vFit; end
% One-parameter fit: v = vMax*(1-(r/(ID/2))^2)
[vMax_fit, parabola_fun] = fit_one_parameter_parabola(rFit(:), vFit(:), ID);
vMean_fit = vMax_fit / 2; % for parabolic profile: vMean = vMax/2

rFine = linspace(0, ID/2, 100);
vFitLine = parabola_fun(vMax_fit, rFine);

fprintf('Parabolic velocity fit: vMax = %.3f cm/s, vMean = %.3f cm/s\n', vMax_fit, vMean_fit);

%% --- Polynomial fit to magnitude vs. radius (blood-only pixels) ---
pMag = polyfit(rGrid(idxBlood).^2, M(idxBlood), 1); % linear in r^2 (parabolic in r)
MFitLine = polyval(pMag, rFine.^2);

%% --- Plot ---
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 38 22]);
hT = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact'); ax = {};

idxSel1 = maskBloodOnly | (maskTissueOnly & rGrid<OD/2);

% velocity vs. radial position
ax{end+1} = nexttile(hT);
plot(rGrid(idxSel1), velPD(idxSel1), '.', 'Color', [0.5 0.5 0.8]);
hold on
plot(rFine, vFitLine, 'r-', 'LineWidth', 2);
xlabel('off-center position [mm]');
ylabel('velocity (cm/s)');
title(sprintf('velocity profile\nvMax=%.2f, vMean=%.2f cm/s', vMax_fit, vMean_fit));
grid on
legend('data','parabolic fit','Location','north');

% magnitude vs. radial position
ax{end+1} = nexttile(hT);
plot(rGrid(idxSel1), M(idxSel1), '.', 'Color', [0.5 0.8 0.5]);
hold on
plot(rFine, MFitLine, 'r-', 'LineWidth', 2);
xlabel('off-center position [mm]');
ylabel('MR signal magnitude [a.u.]');
title('magnitude profile');
grid on
legend('data','polynomial fit','Location','north');

% phase vs. magnitude scatter
ax{end+1} = nexttile(hT);
plot(PD(idxSel1), M(idxSel1), '.', 'Color', [0.8 0.5 0.5]);
xlabel('MR phase difference [rad]');
ylabel('MR signal magnitude [a.u.]');
title('phase vs. magnitude');
grid on

% velocity 2D map
ax{end+1} = nexttile(hT);
vLimPD = max(abs(velPD(:)));
imagesc(PEpos, FEpos, velPD, [-vLimPD vLimPD]); axis image;
ax{end}.Colormap = redblue;
hold on
theta = linspace(0,2*pi,360);
plot(ID/2*cos(theta), ID/2*sin(theta), 'w--', 'LineWidth', 1);
xlabel('PE [mm]'); ylabel('FE [mm]');
title(sprintf('velocity map (venc=%gcm/s)', bestVenc));
ylabel(colorbar, 'velocity [cm/s]');
set(ax{end},'XTick',[],'YTick',[]);

% magnitude 2D map
ax{end+1} = nexttile(hT);
imagesc(PEpos, FEpos, M, [0 max(M(:))]); axis image;
ax{end}.Colormap = gray;
hold on
plot(ID/2*cos(theta), ID/2*sin(theta), 'w--', 'LineWidth', 1);
xlabel('PE [mm]'); ylabel('FE [mm]');
title('magnitude map (venc=inf)');
set(ax{end},'XTick',[],'YTick',[]);

% blood-only mask
ax{end+1} = nexttile(hT);
imagesc(PEpos, FEpos, maskBloodOnly, [0 1]); axis image;
ax{end}.Colormap = gray;
hold on
plot(ID/2*cos(theta), ID/2*sin(theta), 'm', 'LineWidth', 1.5);
title('blood-only mask');
set(ax{end},'XTick',[],'YTick',[]);

drawnow;
if saveThis
    saveas(        f, fullfile(figDir, 'radialProfiles.fig'));
    exportgraphics(f, fullfile(figDir, 'radialProfiles.png'));
    exportgraphics(f, fullfile(figDir, 'radialProfiles.svg'));
    disp('Saved radialProfiles');
end

fprintf('Fit summary for matched simulation:\n  vMax = %.3f cm/s\n  vMean = %.3f cm/s\n  ID = %.2f mm (fixed)\n', vMax_fit, vMean_fit, ID);
disp('dev_slide09_radial done');
