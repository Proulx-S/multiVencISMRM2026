% dev_slide10_11.m -- dev script for in vivo blocks (slides 10-11)
% Run from this file's directory: matlab -batch "run('dev_slide10_11.m')"
% from /scratch/bass/projects/multiVencISMRM2026/tmp/
% Once verified, code gets merged into doIt.m blocks 5-6.

addpath(genpath('/scratch/bass/tools/util'));

inVivoScratch  = '/scratch/bass/projects/multiVencInVivo/tmp';
projectStorage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
figDir         = fullfile(projectStorage, 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

saveThis = 1;

%% ROI coordinates (from multiVencInVivo/doIt.m lines 378-424)
subRoiList = {};
% sub-01: 6 vessels
subRoiList{end+1} = struct();
subRoiList{end}(1).roiY = [37 47];   subRoiList{end}(1).roiX = [87 92];
subRoiList{end}(2).roiY = [158 164]; subRoiList{end}(2).roiX = [90 92];
subRoiList{end}(3).roiY = [91 94];   subRoiList{end}(3).roiX = [88 91];
subRoiList{end}(4).roiY = [103 107]; subRoiList{end}(4).roiX = [77 81];
subRoiList{end}(5).roiY = [100 103]; subRoiList{end}(5).roiX = [139 141];
subRoiList{end}(6).roiY = [130 134]; subRoiList{end}(6).roiX = [50 53];
% sub-02: 5 vessels
subRoiList{end+1} = struct();
subRoiList{end}(1).roiY = [235 242]; subRoiList{end}(1).roiX = [140 144];
subRoiList{end}(2).roiY = [228 232]; subRoiList{end}(2).roiX = [55 59];
subRoiList{end}(3).roiY = [199 203]; subRoiList{end}(3).roiX = [136 140];
subRoiList{end}(4).roiY = [216 219]; subRoiList{end}(4).roiX = [52 54];
subRoiList{end}(5).roiY = [163 169]; subRoiList{end}(5).roiX = [91 94];

subNames = {'sub-01','sub-02'};

for s = 1:2
    subFile = fullfile(inVivoScratch, [subNames{s} '.mat']);
    load(subFile,'img','imgInfo','refImgAv');
    disp(['Loaded ' subNames{s}]);

    %% ROI overlay figure
    hF = figure('MenuBar','none','ToolBar','none','Units','normalized','Position',[0 0 1 1]);
    hIm = imagesc(mean(abs(refImgAv),7));
    for roiIdx = 1:length(subRoiList{s})
        roiY = subRoiList{s}(roiIdx).roiY;
        roiX = subRoiList{s}(roiIdx).roiX;
        hold on
        hBox(roiIdx)  = plot([roiX(1)-0.5 roiX(2)+0.5 roiX(2)+0.5 roiX(1)-0.5 roiX(1)-0.5], ...
                             [roiY(1)-0.5 roiY(1)-0.5 roiY(2)+0.5 roiY(2)+0.5 roiY(1)-0.5], 'c');
        hText(roiIdx) = text(roiX(1),roiY(1),sprintf('roi%d',roiIdx),'Color','r','FontSize',12,'FontWeight','bold');
    end
    ax = gca; axis image;
    set(ax,'XTick',[],'YTick',[],'Colormap',gray,'DataAspectRatio',[imgInfo.res 1]);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.CLim = [0, 1/2*max(hIm.CData(:))];
    set(hBox,'LineWidth',0.5);
    drawnow;
    if saveThis
        exportgraphics(hF, fullfile(figDir, [subNames{s} '-roiOverlay.png']));
        exportgraphics(hF, fullfile(figDir, [subNames{s} '-roiOverlay.svg']));
        disp(['  saved roiOverlay for ' subNames{s}]);
    end
    close(hF);

    %% Per-vessel spiral figures
    for roiIdx = 1:length(subRoiList{s})
        roiY = subRoiList{s}(roiIdx).roiY;
        roiX = subRoiList{s}(roiIdx).roiX;

        % Extract vessel ROI
        trj = img(roiY(1):roiY(2), roiX(1):roiX(2), :,:,:,:,:,:,:,:,:,:,:,:,:,:);

        % Remove reference phase per run
        runIdxList = unique(imgInfo.runIdx);
        for runIdx = 1:length(runIdxList)
            idx = squeeze(imgInfo.runIdx==runIdxList(runIdx) & imgInfo.vencList==inf);
            refPhase = angle(mean(trj(:,:,:,:,:,:,idx,:,:,:,:,:,:,:,:,:), [7 11]));
            idx2 = squeeze(imgInfo.runIdx==runIdxList(runIdx));
            trj(:,:,:,:,:,:,idx2,:,:,:,:,:,:,:,:,:) = ...
                trj(:,:,:,:,:,:,idx2,:,:,:,:,:,:,:,:,:) ./ exp(1i*refPhase);
        end

        % Average over spatial dims and cardiac, reshape to [venc*run, time]
        trj     = permute(mean(trj,[1 2]),[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
        trj     = trj ./ abs(mean(trj(1:2,:),[1 2]));
        trjVenc = permute(imgInfo.vencList,[7 11 1 2 3 4 5 6 8 9 10 12 13 14 15 16]);
        trjVencLabel = replace(cellstr(num2str(trjVenc)),' ','');

        % Credible mean per venc
        d = {};
        for vencIdx = 1:size(trj,1)
            d{vencIdx,1} = [real(trj(vencIdx,:)); imag(trj(vencIdx,:))]';
        end
        prob = credibleMean2d(d);

        % Plot
        hF2 = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 20 20]);
        cMap_venc = sort(unique(trjVenc),'descend');
        cMap = jet(length(cMap_venc));
        for dIdx = 1:length(d)
            cMap_idx = cMap_venc==trjVenc(dIdx);
            hPcont(dIdx) = plot(prob{dIdx}.CIcontour,'FaceColor',cMap(cMap_idx,:),'EdgeColor','k');
            hold on
        end
        ax2 = gca;
        set(hPcont,'FaceAlpha',1);
        title(['vessel' num2str(roiIdx,'%02d')]);
        axis image tight
        xLim = xlim; if xLim(1)>0; xLim(1)=0; end; if xLim(2)<0; xLim(2)=0; end;
        yLim = ylim; if yLim(1)>0; yLim(1)=0; end; if yLim(2)<0; yLim(2)=0; end;
        dLim = max(diff(xLim),diff(yLim)).*0.03;
        xLim = xLim + [-dLim dLim];
        yLim = yLim + [-dLim dLim];
        set(ax2,'XLim',xLim,'YLim',yLim);
        uistack(xline(0,'w'),'bottom');
        uistack(yline(0,'w'),'bottom');
        grid on
        xlabel('Re','FontName','Times New Roman');
        ylabel('Im','FontName','Times New Roman');
        ax2.Color = 'k';
        dTick = min(mean(diff(ax2.XTick)),mean(diff(ax2.YTick)));
        xTick = round(xLim./dTick).*dTick;
        xTick = linspace(xTick(1),xTick(2),range(xTick)/dTick+1);
        yTick = round(yLim./dTick).*dTick;
        yTick = linspace(yTick(1),yTick(2),range(yTick)/dTick+1);
        set(ax2,'XTick',xTick,'YTick',yTick);
        theta = linspace(0,2*pi,100);
        x = abs(mean(trj(trjVenc==inf,:),[1 2]))*cos(theta);
        y = abs(mean(trj(trjVenc==inf,:),[1 2]))*sin(theta);
        uistack(plot(x,y,'w'),'bottom');
        title(legend(hPcont,trjVencLabel),'V_{enc} (cm/s)');
        drawnow;
        if saveThis
            exportgraphics(hF2, fullfile(figDir, [subNames{s} '_vessel-' num2str(roiIdx,'%02d') '.png']));
            exportgraphics(hF2, fullfile(figDir, [subNames{s} '_vessel-' num2str(roiIdx,'%02d') '.svg']));
            disp(['  saved vessel ' num2str(roiIdx,'%02d') ' for ' subNames{s}]);
        end
        close(hF2);
        clear hPcont hBox hText
    end
end
disp('dev_slide10_11 done');
