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
tool = 'util'; repoURL = 'https://github.com/Proulx-S/util.git'; subTool = ''; branch = 'dev-inflow';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
tool = 'pcMRAsim'; repoURL = 'https://github.com/Proulx-S/pcMRAsim.git'; subTool = ''; branch = 'dev-inflow';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
tool = 'multiVencSim'; repoURL = 'https://github.com/Proulx-S/multiVencSim.git'; subTool = ''; branch = 'dev-inflow';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
% pcMRAsim must be re-added last: multiVencSim carries its own old runSim.m which shadows pcMRAsim's
addpath(genpath(fullfile(toolDir, 'pcMRAsim')));
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

% Compute coordinates around center of mass
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, size(M,1)*FEspacing-FEspacing/2, size(M,1));
PEpos = linspace(PEspacing/2, size(M,2)*PEspacing-PEspacing/2, size(M,2));
[FEgrid, PEgrid] = ndgrid(FEpos, PEpos);
total = sum(M(:));
com(1) = sum(FEgrid(:) .* M(:)) / total;
com(2) = sum(PEgrid(:) .* M(:)) / total;
FEgrid = FEgrid - com(1);
FEpos  = FEpos  - com(1);
PEgrid = PEgrid - com(2);
PEpos  = PEpos  - com(2);
rGrid = sqrt(PEgrid.^2+FEgrid.^2);
pGrid = -atan2(FEgrid, PEgrid);  % theta=0 → +PE ("right" in imagesc display)

theta = linspace(0, 2*pi, 360);

clear dFE dPE d_far d_near M com total
%% %%%%%%%%%%%%%%%%%


% [Sections 1–10 from doIt.m: poster/background figures — not reproduced here.
%  Convention: each section lives in `if 0 ... end`, opens with `saveThis = 0;`,
%  has its own `sec<N>fig` folder, and ends with `%% [delimiter]` + `end % section N`.
%  Activate a section by changing `if 0` → `if 1`.]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in vivo data -- sub-01 and sub-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI coordinates (vessels defined in multiVencInVivo project)
inVivoSubRoiList = {};
% sub-01: 6 vessels  (bestVenc from multiVencInVivo/doIt.m per-vessel selection)
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY = [37 47];   inVivoSubRoiList{end}(1).roiX = [87 92];  inVivoSubRoiList{end}(1).bestVenc = 10;
inVivoSubRoiList{end}(2).roiY = [158 164]; inVivoSubRoiList{end}(2).roiX = [90 92];  inVivoSubRoiList{end}(2).bestVenc = 13;
inVivoSubRoiList{end}(3).roiY = [91 94];   inVivoSubRoiList{end}(3).roiX = [88 91];  inVivoSubRoiList{end}(3).bestVenc = 20;
inVivoSubRoiList{end}(4).roiY = [103 107]; inVivoSubRoiList{end}(4).roiX = [77 81];  inVivoSubRoiList{end}(4).bestVenc = 20;
inVivoSubRoiList{end}(5).roiY = [100 103]; inVivoSubRoiList{end}(5).roiX = [139 141]; inVivoSubRoiList{end}(5).bestVenc = 5;
inVivoSubRoiList{end}(6).roiY = [130 134]; inVivoSubRoiList{end}(6).roiX = [50 53];  inVivoSubRoiList{end}(6).bestVenc = 5;
% sub-02: 5 vessels
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY = [235 242]; inVivoSubRoiList{end}(1).roiX = [140 144]; inVivoSubRoiList{end}(1).bestVenc = 8;
inVivoSubRoiList{end}(2).roiY = [228 232]; inVivoSubRoiList{end}(2).roiX = [55 59];   inVivoSubRoiList{end}(2).bestVenc = 40;
inVivoSubRoiList{end}(3).roiY = [199 203]; inVivoSubRoiList{end}(3).roiX = [136 140]; inVivoSubRoiList{end}(3).bestVenc = 7;
inVivoSubRoiList{end}(4).roiY = [216 219]; inVivoSubRoiList{end}(4).roiX = [52 54];   inVivoSubRoiList{end}(4).bestVenc = 4;
inVivoSubRoiList{end}(5).roiY = [163 169]; inVivoSubRoiList{end}(5).roiX = [91 94];   inVivoSubRoiList{end}(5).bestVenc = 5;

inVivoSubNames   = {'sub-01','sub-02'};
inVivoScratch    = fullfile(fileparts(info.project.code), 'multiVencInVivo', 'tmp');
inVivoSubData    = cell(1,2);
for s = 1:2
    subFile = fullfile(inVivoScratch, [inVivoSubNames{s} '.mat']);
    inVivoSubData{s} = load(subFile, 'img', 'imgInfo', 'refImgAv');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% return


if 0
saveThis = 0;
%% 11 - [placeholder: Fit A — see doIt.m sec11 for full implementation]
%% %%%
end % section 11


if 0
saveThis = 0;
%% 12 - [placeholder: Fit B — see doIt.m sec12 for full implementation]
%% %%%
end % section 12


if 1
saveThis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13 - inflow model (phantom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec13fig = fullfile(info.project.figures, '13-inflow-model');
if ~exist(sec13fig,'dir'); mkdir(sec13fig); end

% --- MRI parameters (phantom 03: BEAT-FQ sequence, 3T) ---
p_ph_def = runSim;
pMri_ph  = p_ph_def.pMri;
pMri_ph.fieldStrength  = 3;
pMri_ph.species        = 'phantom';
pMri_ph.sliceThickness = 2.2;          % mm
pMri_ph.TR             = 75.90/(5+1)/1000;  % s  (75.90 ms / 6 echoes)
pMri_ph.TE             = 9.8/1000;     % s
pMri_ph.FA             = 50;           % deg
pMri_ph.venc.method    = 'FVEmono';
pMri_ph.venc.FVEbw     = 100;          % cm/s
p_ph = runSim(p_ph_def.pVessel, p_ph_def.pSim, pMri_ph);
pMri_ph   = p_ph.pMri;
pRelax_ph = pMri_ph.relax.blood;

% --- Extract phantom data ---
cFlow_ph   = squeeze(mean(data(:,:,dataVenc==inf),    3));
cBest_ph   = squeeze(mean(data(:,:,dataVenc==bestVenc),3));
cNoFlow_ph = squeeze(mean(dataNoFlow(:,:,dataVenc==inf),3));
M_ph       = abs(cFlow_ph);
vFlow_ph   = phase2vel(angle(cBest_ph), vencToM1(bestVenc));

maskBlood_ph = M_ph > 0.30 * max(M_ph(:));
r_blood = rGrid(maskBlood_ph);   p_blood = pGrid(maskBlood_ph);
m_blood = double(M_ph(maskBlood_ph));    v_blood = double(vFlow_ph(maskBlood_ph));
m_noflow_blood = double(abs(cNoFlow_ph(maskBlood_ph)));
R_ph = ID/2;   % inner lumen radius [mm]

% --- Physics reference at v=0 (fully saturated) ---
Mz_v0  = getMz_ss(pMri_ph, pRelax_ph, 0);
Mxy_v0 = double(getMxy_ss(Mz_v0, pMri_ph, pRelax_ph));
A_init = mean(m_noflow_blood) / Mxy_v0;

% --- Pre-realize complex noise (fixed before optimization; same realization every iteration) ---
rng(0);
noise_grid   = randn(size(M_ph)) + 1i * randn(size(M_ph));   % for flow data
noise_noflow = randn(size(M_ph)) + 1i * randn(size(M_ph));   % for noFlow (independent)

% --- Initial estimates for new parameters ---
WT_init          = double(OD/2 - ID/2);                          % physical tube wall [mm]
mask_wall_init   = rGrid > R_ph & rGrid <= OD/2;
mask_tissue_init = rGrid >  OD/2;
sigma_n_init  = double(mean(M_ph(mask_wall_init))) * sqrt(2/pi); % Rayleigh E → Gaussian sigma
S_tissue_init = double(mean(M_ph(mask_tissue_init)));
if isnan(S_tissue_init) || isempty(S_tissue_init)
    S_tissue_init = double(mean(M_ph(:))) * 0.3;
end

% --- Fit bounds and initial values ---
% theta = [Vmax, R, tx, ty, A, FEoffset, PEoffset, WT, S_tissue, sigma_n]
sv        = std(v_blood);
sm        = double(std(M_ph(:)));   % all pixels for three-compartment fit
Vmax_init = max(abs(v_blood));
lb_13 = double([0,    1e-6, -0.7, -0.7, 0,    -PEspacing,  -FEspacing,  0,      0,              0  ]);
ub_13 = double([inf,  ID,    0.7,  0.7, inf,   PEspacing,   FEspacing,   OD,    inf,             inf]);
theta0_13 = double([Vmax_init, R_ph, 0, 0, A_init, 0, 0, WT_init, S_tissue_init, sigma_n_init]);
theta0_13 = min(max(theta0_13, lb_13), ub_13);
opts13 = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',3e4,'FunctionTolerance',1e-9);

% --- Fit 13a: three-compartment + noise, with noFlow ---
f_res_a = @(th) residuals_inflow13_full(th, FEgrid, PEgrid, ...
    double(M_ph), double(vFlow_ph), maskBlood_ph, double(abs(cNoFlow_ph)), ...
    noise_grid, noise_noflow, pMri_ph, pRelax_ph, sv, sm);
theta_13a = lsqnonlin(f_res_a, theta0_13, lb_13, ub_13, opts13);

% --- Fit 13b: three-compartment + noise, without noFlow ---
f_res_b = @(th) residuals_inflow13_full(th, FEgrid, PEgrid, ...
    double(M_ph), double(vFlow_ph), maskBlood_ph, [], ...
    noise_grid, noise_noflow, pMri_ph, pRelax_ph, sv, sm);
theta_13b = lsqnonlin(f_res_b, theta_13a, lb_13, ub_13, opts13);

% --- Derive physical parameters ---
[thetaDeg_a, AR_a, alphaDeg_a, pMri_eff_a] = vessel_angle_params13(theta_13a(3), theta_13a(4), pMri_ph);
[thetaDeg_b, AR_b, alphaDeg_b, pMri_eff_b] = vessel_angle_params13(theta_13b(3), theta_13b(4), pMri_ph);
alpha_b_rad = alphaDeg_b * pi/180;
R_b  = theta_13b(2);
WT_b = theta_13b(8);

% --- Complex domain: data trajectory ---
finiteVencs_ph = sort(unique(dataVenc(~isinf(dataVenc))));
m1_meas_ph = arrayfun(@vencToM1, finiteVencs_ph);

all_vencs_ph = [inf; finiteVencs_ph(:)];
trj_ph = zeros(numel(all_vencs_ph), 1);
trj_ph(1) = mean(cFlow_ph(maskBlood_ph));
for kk = 1:numel(finiteVencs_ph)
    cVenc_kk = squeeze(mean(data(:,:,dataVenc==finiteVencs_ph(kk)), 3));
    trj_ph(1+kk) = mean(cVenc_kk(maskBlood_ph));
end
trj_ph_n = trj_ph / abs(trj_ph(1));

% --- Figure prep ---
rGridOff_ph = sqrt((FEgrid - theta_13b(6)).^2 + (PEgrid - theta_13b(7)).^2);
r_max_plt   = max(rGridOff_ph(:)) * 1.02;
r_plt       = linspace(0, r_max_plt, 300);

% 1D radial profile segments — three compartments (Fit 13b, along p=0)
v1D_b = @(r) velocity_func_ellipse(r, zeros(size(r)), theta_13b(1), R_b, AR_b, alpha_b_rad, 0, 0);
r_lumen_plt  = r_plt(r_plt <  R_b);
r_wall_plt   = r_plt(r_plt >= R_b & r_plt < R_b + WT_b);
r_tissue_plt = r_plt(r_plt >= R_b + WT_b);
m_r_lumen_b  = inflowMag13(v1D_b(r_lumen_plt), theta_13b(5), pMri_eff_b, pRelax_ph);
m_r_wall_b   = theta_13b(10) * sqrt(pi/2) * ones(size(r_wall_plt));  % Rayleigh E[|noise|]
m_r_tissue_b = theta_13b(9)  * ones(size(r_tissue_plt));

% m(v) inflow curve (blood pixels only)
v_plt   = linspace(0, theta_13b(1)*1.1, 200);
m_plt_b = inflowMag13(v_plt, theta_13b(5), pMri_eff_b, pRelax_ph);

% Ellipse overlays — inner and outer wall (Fit 13b)
t_c    = linspace(0, 2*pi, 300);
PE_in  = R_b            .* cos(t_c);   FE_in  = (R_b/AR_b)          .* sin(t_c);
PE_out = (R_b + WT_b)   .* cos(t_c);   FE_out = (R_b/AR_b + WT_b)   .* sin(t_c);
cx_in  = PE_in.*cos(alpha_b_rad)  - FE_in.*sin(alpha_b_rad)  + theta_13b(7);
cy_in  = PE_in.*sin(alpha_b_rad)  + FE_in.*cos(alpha_b_rad)  + theta_13b(6);
cx_out = PE_out.*cos(alpha_b_rad) - FE_out.*sin(alpha_b_rad) + theta_13b(7);
cy_out = PE_out.*sin(alpha_b_rad) + FE_out.*cos(alpha_b_rad) + theta_13b(6);

% Three-compartment predicted magnitudes for scatter plots (all voxels, Fit 13b)
v_pred_full = velocity_func_ellipse(rGrid, pGrid, theta_13b(1), R_b, AR_b, alpha_b_rad, theta_13b(6), theta_13b(7));
dPE_f = PEgrid - theta_13b(7);   dFE_f = FEgrid - theta_13b(6);
r_v_f = sqrt(dPE_f.^2 + dFE_f.^2);
uPE_f = dPE_f./max(r_v_f,eps);   uFE_f = dFE_f./max(r_v_f,eps);
Ae_f  = uPE_f.*cos(alpha_b_rad) + uFE_f.*sin(alpha_b_rad);
Be_f  = -uPE_f.*sin(alpha_b_rad) + uFE_f.*cos(alpha_b_rad);
Reff_in_f  = R_b ./ sqrt(max(Ae_f.^2 + AR_b^2.*Be_f.^2, eps));
lumen_f    = v_pred_full > 0;
wall_f     = ~lumen_f & (r_v_f < Reff_in_f + WT_b);
Mxy_f      = zeros(size(M_ph));
Mxy_f(lumen_f) = inflowMag13(v_pred_full(lumen_f), theta_13b(5), pMri_eff_b, pRelax_ph);
Mxy_f(~lumen_f & ~wall_f) = theta_13b(9);
m_pred_full    = sqrt((Mxy_f + theta_13b(10)*real(noise_grid)).^2 + (theta_13b(10)*imag(noise_grid)).^2);
phase_pred_full = double(pi * v_pred_full / bestVenc);
phase_meas_full = double(angle(cBest_ph));

% Analytical complex trajectory from Fit 13b — blood pixels only, phase-aligned to data
% phase = gamma_Hz * M1[T·s²/m] * v[m/s]
gamma_hz  = 2.6752218708e8 / (2*pi);
v_b_ms    = double(v_pred_full(maskBlood_ph)) / 100;   % cm/s → m/s
Mxy_b_cd  = double(Mxy_f(maskBlood_ph));
m1_traj   = linspace(0, vencToM1(min(finiteVencs_ph)), 600)';   % M1=0 → max M1
phases_cd = gamma_hz .* m1_traj .* v_b_ms(:)';            % 600 × N_blood
I_pred_cd = mean(Mxy_b_cd(:)' .* exp(1j .* phases_cd), 2);% 600 × 1
I_pred_cd_n = I_pred_cd / abs(I_pred_cd(1));
I_pred_cd_n = I_pred_cd_n * exp(1j * angle(trj_ph_n(1))); % align to data background phase

% --- Sub-grid partial-volume fractions per voxel (7×7 sub-points) ---
n_sub = 7;
[dfe_sg, dpe_sg] = ndgrid(linspace(-0.5,0.5,n_sub)*FEspacing, linspace(-0.5,0.5,n_sub)*PEspacing);
dfe_sg = dfe_sg(:)';   dpe_sg = dpe_sg(:)';   % 1 × n_sub²
fe_sg  = FEgrid(:) + dfe_sg;   % nPix × n_sub²
pe_sg  = PEgrid(:) + dpe_sg;
r_sg   = sqrt(fe_sg.^2 + pe_sg.^2);
p_sg   = -atan2(fe_sg, pe_sg);
v_sg   = velocity_func_ellipse(r_sg(:), p_sg(:), theta_13b(1), R_b, AR_b, alpha_b_rad, theta_13b(6), theta_13b(7));
v_sg   = reshape(v_sg, numel(M_ph), n_sub^2);
dPE_sg = pe_sg - theta_13b(7);   dFE_sg = fe_sg - theta_13b(6);
rv_sg  = sqrt(dPE_sg.^2 + dFE_sg.^2);
uPE_sg = dPE_sg ./ max(rv_sg, eps);   uFE_sg = dFE_sg ./ max(rv_sg, eps);
Ae_sg  = uPE_sg.*cos(alpha_b_rad) + uFE_sg.*sin(alpha_b_rad);
Be_sg  = -uPE_sg.*sin(alpha_b_rad) + uFE_sg.*cos(alpha_b_rad);
Reff_sg = reshape(R_b ./ sqrt(max(Ae_sg.^2 + AR_b^2.*Be_sg.^2, eps)), numel(M_ph), n_sub^2);
rv_sg   = reshape(rv_sg, numel(M_ph), n_sub^2);
f_lum = mean(v_sg > 0, 2);
f_wal = mean(v_sg == 0 & rv_sg < Reff_sg + WT_b, 2);
f_tis = max(0, 1 - f_lum - f_wal);
[~, comp_dom] = max([f_lum, f_wal, f_tis], [], 2);
comp_map = reshape(comp_dom, size(M_ph));   % 1=lumen  2=wall  3=tissue

clr_l = [0.15 0.85 1.00];   % cyan-blue  → lumen / blood
clr_w = [1.00 0.45 0.05];   % orange     → wall
clr_t = [0.45 0.90 0.45];   % green      → tissue

mks = 4;   mfc = 'w';  mec = 'k';  mlw = 0.5;

f_13 = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 40 16]);
tl_13 = tiledlayout(f_13, 2, 4, 'TileSpacing','compact','Padding','compact');

% (1,1) mag map — inner wall (solid), outer wall (dashed), center +
nexttile(1); imagesc(PEpos, FEpos, M_ph); axis image; colormap(gca,gray); colorbar; hold on;
plot(cx_in,  cy_in,  'r-',  'LineWidth', 1.5);
plot(cx_out, cy_out, 'r--', 'LineWidth', 1.0);
plot(theta_13b(7), theta_13b(6), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
title('mag | venc=\infty | inflow fit'); set(gca,'XTick',[],'YTick',[]);

% (2,1) vel map — same overlays
nexttile(5); imagesc(PEpos, FEpos, vFlow_ph, [-bestVenc bestVenc]); axis image;
colormap(gca, blueBlackRed); colorbar; hold on;
plot(cx_in,  cy_in,  'r-',  'LineWidth', 1.5);
plot(cx_out, cy_out, 'r--', 'LineWidth', 1.0);
plot(theta_13b(7), theta_13b(6), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
title(['vel | venc=' num2str(bestVenc) ' cm/s']); set(gca,'XTick',[],'YTick',[]);

% (1,2) mag radial profile — ALL voxels, three-compartment fit, no grid
nexttile(2);
cm = comp_map(:);  rx = rGridOff_ph(:);  my = double(M_ph(:));
h2l = plot(rx(cm==1), my(cm==1), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_l, 'MarkerEdgeColor',mec, 'LineWidth',mlw); hold on;
h2w = plot(rx(cm==2), my(cm==2), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_w, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
h2t = plot(rx(cm==3), my(cm==3), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_t, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
xline(R_b,        'r-',  'HandleVisibility','off');
xline(R_b + WT_b, 'r--', 'HandleVisibility','off');
h2_fit = plot(r_lumen_plt,  m_r_lumen_b,  'c-', 'LineWidth', 1.5);
         plot(r_wall_plt,   m_r_wall_b,   'c-', 'LineWidth', 1.5, 'HandleVisibility','off');
         plot(r_tissue_plt, m_r_tissue_b, 'c-', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('r [mm]'); ylabel('mag [a.u.]');
legend([h2l h2w h2t h2_fit], {'blood','wall','tissue','fit 13b'}, 'Location','best', 'TextColor','w', 'Color','k');
title('mag radial profile'); set(gca,'Color','k'); axis square;

% (2,2) vel radial profile — ALL voxels, no grid
nexttile(6);
vy = double(vFlow_ph(:));
h6l = plot(rx(cm==1), vy(cm==1), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_l, 'MarkerEdgeColor',mec, 'LineWidth',mlw); hold on;
h6w = plot(rx(cm==2), vy(cm==2), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_w, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
h6t = plot(rx(cm==3), vy(cm==3), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_t, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
xline(R_b,        'r-',  'HandleVisibility','off');
xline(R_b + WT_b, 'r--', 'HandleVisibility','off');
h6_fit = plot(r_lumen_plt,  v1D_b(r_lumen_plt),      'c-', 'LineWidth', 1.5);
         plot(r_wall_plt,   zeros(size(r_wall_plt)),  'c-', 'LineWidth', 1.5, 'HandleVisibility','off');
         plot(r_tissue_plt, zeros(size(r_tissue_plt)),'c-', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('r [mm]'); ylabel('v [cm/s]');
legend([h6l h6w h6t h6_fit], {'blood','wall','tissue','fit 13b'}, 'Location','best', 'TextColor','w', 'Color','k');
title('vel radial profile'); set(gca,'Color','k'); axis square;

% (1,3) pred vs meas magnitude — same x/y range and ticks
nexttile(3);
xy_lim_m = [0, max(max(double(M_ph(:))), max(m_pred_full(:))) * 1.02];
plot(xy_lim_m, xy_lim_m, '--', 'Color', [.35 .35 .35], 'LineWidth', 0.8); hold on;
mx = double(M_ph(:));  my3 = m_pred_full(:);
plot(mx(cm==1), my3(cm==1), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_l, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
plot(mx(cm==2), my3(cm==2), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_w, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
plot(mx(cm==3), my3(cm==3), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_t, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
xlim(xy_lim_m); ylim(xy_lim_m);
xlabel('measured mag [a.u.]'); ylabel('predicted mag [a.u.]');
title('pred vs meas mag (all vox)'); set(gca,'Color','k'); axis square;

% (2,3) complex domain — sim line only (no extra markers), behind data
ax_cd = nexttile(7);
plotComplexDomain(ax_cd, trj_ph_n, all_vencs_ph, 'full', 'markers');
set(findobj(ax_cd,'Type','line','Marker','o'), ...
    'MarkerSize',mks+2, 'MarkerFaceColor',mfc, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
hold(ax_cd,'on');
h_sl = plot(ax_cd, real(I_pred_cd_n), imag(I_pred_cd_n), 'c-', 'LineWidth', 1.5);
uistack(h_sl, 'bottom');
legend(ax_cd, {'phantom','fit 13b'}, 'Location','southwest', 'TextColor','w', 'Color','k');
title(ax_cd, 'complex domain');

% (1,4) inflow model m(v) — blood pixels, no grid
nexttile(4);
cm_b = comp_map(maskBlood_ph);
h4l = plot(v_blood(cm_b==1), m_blood(cm_b==1), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_l, 'MarkerEdgeColor',mec, 'LineWidth',mlw); hold on;
h4w = plot(v_blood(cm_b==2), m_blood(cm_b==2), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_w, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
h4t = plot(v_blood(cm_b==3), m_blood(cm_b==3), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_t, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
yline(mean(m_noflow_blood), 'y:', 'noFlow ref', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
h4_fit = plot(v_plt, m_plt_b, 'c-', 'LineWidth', 1.5);
xlabel('v [cm/s]'); ylabel('mag [a.u.]');
h4_handles = []; h4_labels = {};
if ~isempty(v_blood(cm_b==1)); h4_handles(end+1)=h4l; h4_labels{end+1}='blood'; end
if ~isempty(v_blood(cm_b==2)); h4_handles(end+1)=h4w; h4_labels{end+1}='wall';  end
if ~isempty(v_blood(cm_b==3)); h4_handles(end+1)=h4t; h4_labels{end+1}='tissue';end
h4_handles(end+1) = h4_fit;  h4_labels{end+1} = 'fit 13b';
legend(h4_handles, h4_labels, 'Location','best', 'TextColor','w', 'Color','k');
title('inflow model m(v)'); set(gca,'Color','k'); axis square;

% (2,4) pred vs meas phase — same x/y range and ticks
nexttile(8);
plot([-pi pi], [-pi pi], '--', 'Color', [.35 .35 .35], 'LineWidth', 0.8); hold on;
pmx = phase_meas_full(:);  ppy = phase_pred_full(:);
plot(pmx(cm==1), ppy(cm==1), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_l, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
plot(pmx(cm==2), ppy(cm==2), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_w, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
plot(pmx(cm==3), ppy(cm==3), 'o', 'LineStyle','none', 'MarkerSize',mks, 'MarkerFaceColor',clr_t, 'MarkerEdgeColor',mec, 'LineWidth',mlw);
xlim([-pi pi]); ylim([-pi pi]);
set(gca, 'XTick', -pi:pi/2:pi, 'YTick', -pi:pi/2:pi);
xlabel('measured phase [rad]'); ylabel('predicted phase [rad]');
title('pred vs meas phase (all vox)'); set(gca,'Color','k'); axis square;

set(findall(f_13,'Type','axes'),'FontSize',12);
set(findall(f_13,'Type','text'),'FontSize',8);
for k_ = 1:8
    text(nexttile(tl_13,k_), 0.02, 0.97, char('A'+k_-1), ...
        'Units','normalized', 'FontSize',13, 'FontWeight','bold', ...
        'Color','w', 'VerticalAlignment','top', 'HorizontalAlignment','left');
end
if saveThis
    exportgraphics(f_13, fullfile(sec13fig,'phantom_inflow.png'));
    hMkr_=findobj(f_13,'Marker','o'); origMEC_=get(hMkr_,{'MarkerEdgeColor'});
    arrayfun(@(h) set(h,'MarkerEdgeColor',h.MarkerFaceColor), hMkr_); set(hMkr_,'Marker','.');
    exportgraphics(f_13, fullfile(sec13fig,'phantom_inflow.svg'));
    set(hMkr_,'Marker','o'); set(hMkr_,{'MarkerEdgeColor'},origMEC_);
end
close(f_13);

% --- fitInfo and parameter reports ---
paramNames = {'Vmax','R','tx','ty','A','FEoffset','PEoffset','WT','S_tissue','sigma_n'};
paramUnits = {'cm/s','mm','','','a.u.','mm','mm','mm','a.u.','a.u.'};
paramFixed = zeros(1,10);

fitInfo_13a.names  = paramNames;   fitInfo_13a.units  = paramUnits;
fitInfo_13a.fixed  = paramFixed;   fitInfo_13a.lb     = lb_13;   fitInfo_13a.ub = ub_13;
fitInfo_13a.theta0 = theta0_13;    fitInfo_13a.theta  = theta_13a;
fitInfo_13a.derived.names  = {'theta_vessel_deg','AR','alpha_deg'};
fitInfo_13a.derived.units  = {'deg','','deg'};
fitInfo_13a.derived.theta0 = [0, 1, 0];
fitInfo_13a.derived.theta  = [thetaDeg_a, AR_a, alphaDeg_a];
fitInfo_13a.init_notes = { ...
    'max|v_blood| from best-VENC image', ...
    'ID/2 (inner lumen radius)', ...
    '0 (perpendicular vessel)', '0 (perpendicular vessel)', ...
    'mean(noFlow_blood)/getMxy_ss(v=0)', '0', '0', ...
    'OD/2-ID/2 (physical tube wall)', ...
    'mean(pixels outside OD)', ...
    'mean(wall pixels)*sqrt(2/pi)'};

fitInfo_13b           = fitInfo_13a;
fitInfo_13b.theta0    = theta_13a;
fitInfo_13b.theta     = theta_13b;
fitInfo_13b.derived.theta = [thetaDeg_b, AR_b, alphaDeg_b];
fitInfo_13b.init_notes = repmat({'Fit 13a final'}, 1, 10);

N_ph = numel(M_ph);
cost13.equation  = ['\mathbf{r} = \begin{bmatrix}' ...
    '(m_i^{\mathrm{meas}} - |\hat{M}_{xy,i} + \sigma_n \eta_i|)/\sigma_m \\' ...
    '(v_i^{\mathrm{meas}} - v_i^{\mathrm{pred}})/\sigma_v\end{bmatrix}'];
cost13.predicted = '|\hat{M}_{xy,i} + \sigma_n \eta_i|,\quad \eta_i \sim \mathcal{CN}(0,1),\text{ pre-realized}';
cost13.normDesc  = ['$\sigma_m = \mathrm{std}(M^{\mathrm{meas}}_{\mathrm{all}})$, ' ...
    '$\sigma_v = \mathrm{std}(v^{\mathrm{meas}}_{\mathrm{blood}})$. ' ...
    'Three compartments: lumen (inflow model), wall ($\hat{M}_{xy}=0$, WT free), tissue ($S_t$ free). ' ...
    'Noise $\eta_i$ pre-realized once; $\sigma_n$ scales it each iteration.'];
cost13.magModel  = ['$$m_{\mathrm{phys}}(v) = |\hat{M}_{xy}(v) + \sigma_n \eta|$$\n\n' ...
    'Lumen: $\hat{M}_{xy} = A \cdot M_{xy}^{\mathrm{ss}}(\tau(v))$. ' ...
    'Wall: $\hat{M}_{xy}=0$ (fixed). Tissue: $\hat{M}_{xy}=S_t$.'];

dsum_13a = sprintf('%d px (all) + %d vel + %d noFlow. Total: %d.', ...
    N_ph, sum(maskBlood_ph(:)), N_ph, 2*N_ph + sum(maskBlood_ph(:)));
dsum_13b = sprintf('%d px (all) + %d vel. Total: %d.', ...
    N_ph, sum(maskBlood_ph(:)), N_ph + sum(maskBlood_ph(:)));

writeFitParamsMd(fullfile(sec13fig,'phantom_fitParams_13a.md'), ...
    '13a — Inflow 3-compartment (phantom, with noFlow)', fitInfo_13a, cost13, dsum_13a);
writeFitParamsMd(fullfile(sec13fig,'phantom_fitParams_13b.md'), ...
    '13b — Inflow 3-compartment (phantom, no noFlow)',   fitInfo_13b, cost13, dsum_13b);

%% %%%%%%%%%%%%%%%%%%
end % section 13


% =========================================================================
% Local functions
% =========================================================================

function res = residuals_inflow13_full(theta, FEgrid, PEgrid, m_meas, v_meas, mask_vel, m_noflow, noise_grid, noise_noflow, pMri_base, pRelax, sv, sm)
% theta = [Vmax, R, tx, ty, A, FEoffset, PEoffset, WT, S_tissue, sigma_n]
% m_noflow: [] → no noFlow; full image → include noFlow magnitude residuals
Vmax=theta(1); R_=theta(2); tx=theta(3); ty=theta(4); A=theta(5);
FEoffset=theta(6); PEoffset=theta(7); WT=theta(8); S_tissue=theta(9); sigma_n=theta(10);

cosT  = sqrt(max(0, 1 - tx^2 - ty^2));
AR    = 1 / max(cosT, 1e-6);
alpha = atan2(ty, tx);
pMri_eff = pMri_base;
pMri_eff.sliceThickness = pMri_base.sliceThickness / max(cosT, 1e-6);

% Velocity for all pixels
r_all = sqrt(FEgrid.^2 + PEgrid.^2);
p_all = -atan2(FEgrid, PEgrid);
v_pred = velocity_func_ellipse(r_all, p_all, Vmax, R_, AR, alpha, FEoffset, PEoffset);

% Compartment assignment
dPE = PEgrid - PEoffset;   dFE = FEgrid - FEoffset;
r_v = sqrt(dPE.^2 + dFE.^2);
uPE = dPE./max(r_v,eps);   uFE = dFE./max(r_v,eps);
Ae  = uPE.*cos(alpha) + uFE.*sin(alpha);
Be  = -uPE.*sin(alpha) + uFE.*cos(alpha);
R_eff_in = R_ ./ sqrt(max(Ae.^2 + AR.^2.*Be.^2, eps));
lumen = v_pred > 0;
wall  = ~lumen & (r_v < R_eff_in + WT);

% Predicted Mxy
Mxy = zeros(size(FEgrid));
Mxy(lumen) = inflowMag13(v_pred(lumen), A, pMri_eff, pRelax);
% wall: Mxy = 0
Mxy(~lumen & ~wall) = S_tissue;   % tissue

% Predicted magnitude with pre-realized noise (sigma_n scales the fixed realization)
m_pred = sqrt((Mxy + sigma_n.*real(noise_grid)).^2 + (sigma_n.*imag(noise_grid)).^2);

% Magnitude residuals — all voxels
res_mag = (m_meas(:) - m_pred(:)) / sm;

% Velocity residuals — blood-masked lumen pixels
mask_v = mask_vel & lumen;
if any(mask_v(:))
    res_vel = (v_meas(mask_v) - v_pred(mask_v)) / sv;
else
    res_vel = [];
end

% noFlow magnitude residuals
if ~isempty(m_noflow)
    Mz_0      = getMz_ss(pMri_eff, pRelax, 0);
    m_nf_lum  = A * double(getMxy_ss(Mz_0, pMri_eff, pRelax));
    Mxy_nf    = zeros(size(FEgrid));
    Mxy_nf(lumen) = m_nf_lum;
    Mxy_nf(~lumen & ~wall) = S_tissue;
    m_pred_nf = sqrt((Mxy_nf + sigma_n.*real(noise_noflow)).^2 + (sigma_n.*imag(noise_noflow)).^2);
    res_nf    = (m_noflow(:) - m_pred_nf(:)) / sm;
else
    res_nf = [];
end

res = double([res_mag; res_vel; res_nf]);
end


% (legacy stub — superseded by residuals_inflow13_full)
function res = residuals_inflow13(theta, r, p, v_meas, m_meas, m_noflow, pMri_base, pRelax, sv, sm)
Vmax=theta(1); R=theta(2); tx=theta(3); ty=theta(4); A=theta(5);
FEoffset=theta(6); PEoffset=theta(7);

cosT  = sqrt(max(0, 1 - tx^2 - ty^2));
AR    = 1 / max(cosT, 1e-6);
alpha = atan2(ty, tx);

pMri_eff = pMri_base;
pMri_eff.sliceThickness = pMri_base.sliceThickness / max(cosT, 1e-6);

v_pred = velocity_func_ellipse(r, p, Vmax, R, AR, alpha, FEoffset, PEoffset);
m_pred = inflowMag13(v_pred, A, pMri_eff, pRelax);
res_v  = (v_meas(:) - v_pred(:)) / sv;
res_m  = (m_meas(:) - m_pred(:)) / sm;

if ~isempty(m_noflow)
    m_nf_pred = inflowMag13(zeros(numel(m_noflow),1), A, pMri_eff, pRelax);
    res_mNF   = (m_noflow(:) - m_nf_pred(:)) / sm;
else
    res_mNF = [];
end
res = double([res_v; res_m; res_mNF]);
end


function m = inflowMag13(v, A, pMri_eff, pRelax)
[Mz, pMri_u] = getMz_ss(pMri_eff, pRelax, double(v(:)'));
m = double(A) * double(getMxy_ss(double(Mz), pMri_u, pRelax));
m = reshape(m, size(v));
end


function [thetaDeg, AR, alphaDeg, pMri_eff] = vessel_angle_params13(tx, ty, pMri_base)
cosT     = sqrt(max(0, 1 - tx^2 - ty^2));
thetaDeg = acosd(cosT);
AR       = 1 / max(cosT, 1e-6);
alphaDeg = atan2(ty, tx) * 180/pi;
pMri_eff = pMri_base;
pMri_eff.sliceThickness = pMri_base.sliceThickness / max(cosT, 1e-6);
end
