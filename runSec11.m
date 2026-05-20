% Run only the environment setup and section 11 (in-vivo profile fits)
clear all; close all; clc;

projectName = 'multiVencISMRM2026';

os   = char(java.lang.System.getProperty('os.name'));
host = char(java.net.InetAddress.getLocalHost.getHostName);
user = char(java.lang.System.getProperty('user.name'));

if strcmp(os,'Linux') && strcmp(host,'takoyaki') && strcmp(user,'sebp')
    storageDrive  = '/local/users/Proulx-S/';
    scratchDrive  = '/scratch/bass/';
    projectCode    = fullfile(scratchDrive, 'projects', projectName);
    projectStorage = fullfile(storageDrive, 'projects', projectName);
    projectScratch = fullfile(scratchDrive, 'projects', projectName, 'tmp');
    toolDir        = fullfile(scratchDrive, 'tools');
else
    mountPoint = '/Users/sebastienproulx/remote/takoyakiLocal';
    storageDrive   = '/Users/sebastienproulx/bass';
    scratchDrive   = '/Users/sebastienproulx/bass';
    projectCode    = fullfile(scratchDrive, 'projects', projectName);
    projectStorage = fullfile(storageDrive, 'projects', projectName);
    projectScratch = fullfile(scratchDrive, 'projects', projectName, 'tmp');
    toolDir        = fullfile(scratchDrive, 'tools');
end

tool = 'util'; repoURL = 'https://github.com/Proulx-S/util.git'; subTool = ''; branch = 'dev-multiVencISMRM2026';
if ~exist(fullfile(toolDir, tool), 'dir'); system(['git clone ' repoURL ' ' fullfile(toolDir, tool)]); end
addpath(genpath(fullfile(toolDir,'util')));
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);
tool = 'pcMRAsim'; repoURL = 'https://github.com/Proulx-S/pcMRAsim.git'; subTool = ''; branch = 'dev-multiVencISMRM2026';
gitClone(repoURL, fullfile(toolDir, tool), subTool, branch);

info.project.code    = projectCode;
info.project.storage = projectStorage;
info.project.scratch = projectScratch;
info.project.figures = fullfile(projectCode, 'figures');

% ROI list (mirrors doIt.m)
inVivoSubRoiList = {};
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY=[37 47];   inVivoSubRoiList{end}(1).roiX=[87 92];  inVivoSubRoiList{end}(1).bestVenc=10;
inVivoSubRoiList{end}(2).roiY=[158 164]; inVivoSubRoiList{end}(2).roiX=[90 92];  inVivoSubRoiList{end}(2).bestVenc=13;
inVivoSubRoiList{end}(3).roiY=[91 94];   inVivoSubRoiList{end}(3).roiX=[88 91];  inVivoSubRoiList{end}(3).bestVenc=20;
inVivoSubRoiList{end}(4).roiY=[103 107]; inVivoSubRoiList{end}(4).roiX=[77 81];  inVivoSubRoiList{end}(4).bestVenc=20;
inVivoSubRoiList{end}(5).roiY=[100 103]; inVivoSubRoiList{end}(5).roiX=[139 141]; inVivoSubRoiList{end}(5).bestVenc=5;
inVivoSubRoiList{end}(6).roiY=[130 134]; inVivoSubRoiList{end}(6).roiX=[50 53];  inVivoSubRoiList{end}(6).bestVenc=5;
inVivoSubRoiList{end+1} = struct();
inVivoSubRoiList{end}(1).roiY=[235 242]; inVivoSubRoiList{end}(1).roiX=[140 144]; inVivoSubRoiList{end}(1).bestVenc=8;
inVivoSubRoiList{end}(2).roiY=[228 232]; inVivoSubRoiList{end}(2).roiX=[55 59];   inVivoSubRoiList{end}(2).bestVenc=40;
inVivoSubRoiList{end}(3).roiY=[199 203]; inVivoSubRoiList{end}(3).roiX=[136 140]; inVivoSubRoiList{end}(3).bestVenc=7;
inVivoSubRoiList{end}(4).roiY=[216 219]; inVivoSubRoiList{end}(4).roiX=[52 54];   inVivoSubRoiList{end}(4).bestVenc=4;
inVivoSubRoiList{end}(5).roiY=[163 169]; inVivoSubRoiList{end}(5).roiX=[91 94];   inVivoSubRoiList{end}(5).bestVenc=5;

inVivoSubNames = {'sub-01','sub-02'};
inVivoScratch  = fullfile(fileparts(info.project.code), 'multiVencInVivo', 'tmp');
inVivoSubData  = cell(1,2);
for s = 1:2
    subFile = fullfile(inVivoScratch, [inVivoSubNames{s} '.mat']);
    inVivoSubData{s} = load(subFile, 'img', 'imgInfo', 'refImgAv');
end

% ---- Run section 11 ----

saveThis = 1;
sec11fig = fullfile(info.project.figures, '11-in-vivo-profile-fits');
if ~exist(sec11fig,'dir'); mkdir(sec11fig); end

p_iv_defaults   = runSim;
pMri_iv         = p_iv_defaults.pMri;
pMri_iv.fieldStrength  = 3;
pMri_iv.species        = 'human';
pMri_iv.sliceThickness = 4;
pMri_iv.TR             = 8e-3;
pMri_iv.FA             = 25;
pMri_iv.TE             = 3e-3;
pMri_iv.venc.method    = 'FVEmono';
pMri_iv.venc.FVEbw     = 100;
p_iv_resolved  = runSim(p_iv_defaults.pVessel, p_iv_defaults.pSim, pMri_iv);
pMri_iv        = p_iv_resolved.pMri;

vBrainRef = 10;
Mz_v0     = getMz_ss(pMri_iv, pMri_iv.relax.blood, 0);
Mxy_v0    = getMxy_ss(Mz_v0,  pMri_iv, pMri_iv.relax.blood);
Mz_vRef   = getMz_ss(pMri_iv, pMri_iv.relax.blood, vBrainRef);
Mxy_vRef  = getMxy_ss(Mz_vRef, pMri_iv, pMri_iv.relax.blood);
B_init_iv = double(Mxy_v0 / Mxy_vRef);
fprintf('B_init_iv = %.3f (Mxy ratio v0/v%d)\n', B_init_iv, vBrainRef);

for s = 1:2
    img      = inVivoSubData{s}.img;
    imgInfo  = inVivoSubData{s}.imgInfo;
    FEspacing_iv = imgInfo.res(1);
    PEspacing_iv = imgInfo.res(2);

    subFigDir = fullfile(sec11fig, inVivoSubNames{s});
    if ~exist(subFigDir,'dir'); mkdir(subFigDir); end

    for roiIdx = 1:length(inVivoSubRoiList{s})
        roiY        = inVivoSubRoiList{s}(roiIdx).roiY;
        roiX        = inVivoSubRoiList{s}(roiIdx).roiX;
        bestVenc_iv = inVivoSubRoiList{s}(roiIdx).bestVenc;
        figName_iv  = [inVivoSubNames{s} sprintf('_vessel-%02d', roiIdx)];
        fprintf('Processing %s...\n', figName_iv);

        roiImg = img(roiY(1):roiY(2), roiX(1):roiX(2), :,:,:,:,:,:,:,:,:,:,:,:,:,:);
        runIdxList_iv = unique(imgInfo.runIdx);
        for rr = 1:length(runIdxList_iv)
            idx_ref = squeeze(imgInfo.runIdx==runIdxList_iv(rr) & imgInfo.vencList==inf);
            refPhase_rr = angle(mean(roiImg(:,:,:,:,:,:,idx_ref,:,:,:,:,:,:,:,:,:), [7 11]));
            idx_run = squeeze(imgInfo.runIdx==runIdxList_iv(rr));
            roiImg(:,:,:,:,:,:,idx_run,:,:,:,:,:,:,:,:,:) = ...
                roiImg(:,:,:,:,:,:,idx_run,:,:,:,:,:,:,:,:,:) ./ exp(1i*refPhase_rr);
        end

        cFlow_iv  = squeeze(mean(roiImg(:,:,:,:,:,:,imgInfo.vencList==inf,        :,:,:,:,:,:,:,:,:),[7 11]));
        cBest_iv  = squeeze(mean(roiImg(:,:,:,:,:,:,imgInfo.vencList==bestVenc_iv,:,:,:,:,:,:,:,:,:),[7 11]));
        vFlow_iv  = phase2vel(angle(cBest_iv), vencToM1(bestVenc_iv));

        nFE_iv = size(cFlow_iv,1);  nPE_iv = size(cFlow_iv,2);
        FEpos_iv = (0:nFE_iv-1) * FEspacing_iv;
        PEpos_iv = (0:nPE_iv-1) * PEspacing_iv;
        [FEgrid_iv, PEgrid_iv] = ndgrid(FEpos_iv, PEpos_iv);
        M_iv   = abs(cFlow_iv);
        tot_iv = sum(M_iv(:));
        com_iv = [sum(FEgrid_iv(:).*M_iv(:))/tot_iv, sum(PEgrid_iv(:).*M_iv(:))/tot_iv];
        FEgrid_iv = FEgrid_iv - com_iv(1);
        PEgrid_iv = PEgrid_iv - com_iv(2);
        FEpos_iv  = FEpos_iv  - com_iv(1);
        PEpos_iv  = PEpos_iv  - com_iv(2);
        rGrid_iv  = sqrt(FEgrid_iv.^2 + PEgrid_iv.^2);
        pGrid_iv  = -atan2(FEgrid_iv, PEgrid_iv);

        maskBlood_iv = M_iv > 0.30 * max(M_iv(:));
        R_iv = min(diff(roiY)*FEspacing_iv, diff(roiX)*PEspacing_iv) / 2;

        [velFit_iv, magFit_iv, velFit1D_iv] = fitMagVelProfile( ...
            rGrid_iv(maskBlood_iv), pGrid_iv(maskBlood_iv), ...
            vFlow_iv(maskBlood_iv), M_iv(maskBlood_iv), ...
            [], R_iv, 'joint', 2, true, B_init_iv);

        fprintf('  Vmax=%.1f cm/s, R=%.2f mm, FEoff=%.2f, PEoff=%.2f\n', ...
            velFit_iv.Vmax, velFit_iv.R, velFit_iv.FEoffset, velFit_iv.PEoffset);

        rGridOff_iv = sqrt((FEgrid_iv - velFit_iv.FEoffset).^2 + ...
                           (PEgrid_iv - velFit_iv.PEoffset).^2);

        % --- Profiles figure ---
        f_prof = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 24 16]);
        tiledlayout(f_prof, 2, 3, 'TileSpacing','compact','Padding','compact');
        ax_pr = {};

        ax_pr{end+1} = nexttile(1);
        imagesc(PEpos_iv, FEpos_iv, M_iv); axis image;
        ax_pr{end}.Colormap = gray; colorbar;
        title(['mag | venc=\infty | ' inVivoSubNames{s} sprintf(' v%02d',roiIdx)]);
        set(ax_pr{end},'XTick',[],'YTick',[]);

        ax_pr{end+1} = nexttile(4);
        imagesc(PEpos_iv, FEpos_iv, vFlow_iv, [-bestVenc_iv bestVenc_iv]); axis image;
        ax_pr{end}.Colormap = blueBlackRed; colorbar;
        title(['vel | venc=' num2str(bestVenc_iv) ' cm/s']);
        set(ax_pr{end},'XTick',[],'YTick',[]);

        r_plt = linspace(0, velFit_iv.R * 1.1, 120);

        ax_pr{end+1} = nexttile(2);
        plot(ax_pr{end}, rGridOff_iv(maskBlood_iv), M_iv(maskBlood_iv), '.', 'Color',[0.6 0.6 0.6]);
        hold(ax_pr{end},'on');
        plot(ax_pr{end}, r_plt, magFit_iv(velFit1D_iv(r_plt)), 'g-', 'LineWidth',1.5);
        xline(ax_pr{end}, velFit_iv.R, 'w--');
        xlabel('r [mm]'); ylabel('mag [a.u.]');
        title('mag radial profile + fit');
        set(ax_pr{end},'Color','k'); grid(ax_pr{end},'on'); axis square;

        ax_pr{end+1} = nexttile(5);
        plot(ax_pr{end}, rGridOff_iv(maskBlood_iv), vFlow_iv(maskBlood_iv), '.', 'Color',[0.6 0.6 0.6]);
        hold(ax_pr{end},'on');
        plot(ax_pr{end}, r_plt, velFit1D_iv(r_plt), 'g-', 'LineWidth',1.5);
        xline(ax_pr{end}, velFit_iv.R, 'w--');
        xlabel('r [mm]'); ylabel('v [cm/s]');
        title('vel radial profile + fit');
        set(ax_pr{end},'Color','k'); grid(ax_pr{end},'on'); axis square;

        ax_pr{end+1} = nexttile(3);
        v_plt = linspace(0, max(abs(vFlow_iv(maskBlood_iv))), 120);
        plot(ax_pr{end}, vFlow_iv(maskBlood_iv), M_iv(maskBlood_iv), '.', 'Color',[0.6 0.6 0.6]);
        hold(ax_pr{end},'on');
        plot(ax_pr{end}, v_plt, magFit_iv(v_plt), 'g-', 'LineWidth',1.5);
        xline(ax_pr{end}, 0, 'w:');
        xlabel('v [cm/s]'); ylabel('mag [a.u.]');
        title('inflow function m(v)');
        set(ax_pr{end},'Color','k'); grid(ax_pr{end},'on'); axis square;

        set(findall(f_prof,'Type','axes'),'FontSize',12);
        set(findall(f_prof,'Type','text'),'FontSize',8);
        if saveThis
            exportgraphics(f_prof, fullfile(subFigDir, [figName_iv '_profiles.png']));
            hMkr_=findobj(f_prof,'Marker','o'); origMEC_=get(hMkr_,{'MarkerEdgeColor'}); arrayfun(@(h) set(h,'MarkerEdgeColor',h.MarkerFaceColor), hMkr_); set(hMkr_,'Marker','.');
            exportgraphics(f_prof, fullfile(subFigDir, [figName_iv '_profiles.svg']));
            set(hMkr_,'Marker','o'); set(hMkr_,{'MarkerEdgeColor'},origMEC_);
        end
        close(f_prof);

        % --- Matched simulation ---
        pSim_iv    = p_iv_defaults.pSim;
        pVessel_iv = p_iv_defaults.pVessel;
        pMri_iv_s  = pMri_iv;

        pSim_iv.voxGrid.fovFE = nFE_iv * FEspacing_iv;
        pSim_iv.voxGrid.fovPE = nPE_iv * PEspacing_iv;
        pSim_iv.voxGrid.matFE = nFE_iv;
        pSim_iv.voxGrid.matPE = nPE_iv;
        pSim_iv.nSpin         = (2^10)^2;
        pSim_iv.gridMode      = 'pseudoVoxel';

        pVessel_iv.ID      = velFit_iv.R * 2;
        pVessel_iv.WT      = FEspacing_iv;
        pVessel_iv.vMean   = velFit_iv.Vmax / 2;
        pVessel_iv.posFE   = velFit_iv.FEoffset;
        pVessel_iv.posPE   = velFit_iv.PEoffset;
        pVessel_iv.profile = 'parabolic1';

        p_tmp       = runSim(pVessel_iv, pSim_iv, pMri_iv_s);
        pSim_iv     = p_tmp.pSim;
        pVessel_iv  = p_tmp.pVessel;
        pMri_iv_s   = p_tmp.pMri;

        [pSim_iv.spinGrid.feGrid, pSim_iv.spinGrid.peGrid] = ...
            ndgrid(pSim_iv.spinGrid.coorFE, pSim_iv.spinGrid.coorPE);
        pSim_iv.spinGrid.rGrid = sqrt(pSim_iv.spinGrid.peGrid.^2 + pSim_iv.spinGrid.feGrid.^2);
        pSim_iv.spinGrid.pGrid = -atan2(pSim_iv.spinGrid.feGrid, pSim_iv.spinGrid.peGrid);

        pVessel_iv.S.lumen = max(0, magFit_iv(velFit_iv( ...
            pSim_iv.spinGrid.rGrid(pVessel_iv.mask.lumen), ...
            pSim_iv.spinGrid.pGrid(pVessel_iv.mask.lumen))));
        pVessel_iv.S.wall     = 0;
        pVessel_iv.S.surround = 0;

        resSim_iv = runSim(pVessel_iv, pSim_iv, pMri_iv_s, [], false);

        % --- Complex domain overlay ---
        trjIV_s   = permute(mean(roiImg,[1 2]),[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
        trjIV_s   = trjIV_s ./ abs(mean(trjIV_s(1:2,:),[1 2]));
        trjVenc_s = permute(imgInfo.vencList,[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);

        I_sim = squeeze(resSim_iv.I);
        I_sim_norm = I_sim / max(abs(I_sim));

        f_cd = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 16 16]);
        ax_cd = axes(f_cd);
        plotComplexDomain(ax_cd, trjIV_s(:), trjVenc_s(:), 'full', 'markers');
        set(findobj(ax_cd,'Type','line','Marker','o'),'MarkerFaceColor','g','MarkerEdgeColor','k');
        hold(ax_cd,'on');
        plot(ax_cd, real(I_sim_norm), imag(I_sim_norm), 'm-', 'LineWidth',1.5);
        legend(ax_cd, {'in vivo','simulation'}, 'Location','best','TextColor','w','Color','k');
        title(ax_cd, [inVivoSubNames{s} sprintf(' vessel-%02d — complex domain', roiIdx)]);
        set(findall(f_cd,'Type','axes'),'FontSize',12);
        set(findall(f_cd,'Type','text'),'FontSize',8);
        if saveThis
            exportgraphics(f_cd, fullfile(subFigDir, [figName_iv '_complexDomain.png']));
            hMkr_=findobj(f_cd,'Marker','o'); origMEC_=get(hMkr_,{'MarkerEdgeColor'}); arrayfun(@(h) set(h,'MarkerEdgeColor',h.MarkerFaceColor), hMkr_); set(hMkr_,'Marker','.');
            exportgraphics(f_cd, fullfile(subFigDir, [figName_iv '_complexDomain.svg']));
            set(hMkr_,'Marker','o'); set(hMkr_,{'MarkerEdgeColor'},origMEC_);
        end
        close(f_cd);
        fprintf('  saved %s\n', figName_iv);
    end
end
disp('Section 11 complete.');
