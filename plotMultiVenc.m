function [ax,hP] = plotMultiVenc(ax,I,venc,run,windowType,cMapName)
if ~exist('ax','var') || isempty(ax); ax = axes; end
if ~exist('windowType','var') || isempty(windowType); windowType = 'full'; end % 'tight' or 'full'
if ~exist('cMapName','var') || isempty(cMapName); cMapName = 'jet'; end

vencList = sort(unique(venc),'descend');
Mnorm = abs(mean(I(venc==inf)));
for vencIdx = 1:size(vencList,1)
    d{vencIdx,1}     = [real(I(venc==vencList(vencIdx))) imag(I(venc==vencList(vencIdx)))] ./ Mnorm;
    dVenc{vencIdx,1} = venc(venc==vencList(vencIdx));
    dVencLabel(vencIdx,1) = unique(dVenc{vencIdx,1});
    dRun{vencIdx,1}  = run( venc==vencList(vencIdx));
end



% cMap = jet;
cMap = eval(cMapName);

M1 = vencToM1(dVencLabel);
cMapM1 = M1-min(M1);
cMapM1 = cMapM1./max(cMapM1);
% cMap = interp1(linspace(0,1,size(cMap,1))',cMap,linspace(0,1,length(d))');
cMap = interp1(linspace(0,1,size(cMap,1))',cMap,cMapM1);

% dProb = credibleMean2d(d);
% for dIdx = 1:length(d)
%     hP(dIdx) = plot(ax,dProb{dIdx}.CIcontour,'FaceColor',cMap(dIdx,:),'EdgeColor','none','FaceAlpha',1); hold on
% end
for dIdx = 1:length(d)
    hP(dIdx) = plot(ax,mean(d{dIdx,1}(:,1),1),mean(d{dIdx,1}(:,2),1),'o','MarkerFaceColor',cMap(dIdx,:),'MarkerEdgeColor','w'); hold on
end
set(hP,'MarkerSize',10,'LineWidth',1);
axis image tight
switch windowType
    case 'full'
        xLim = [-1 1].*max(abs([xlim ylim]));
        yLim = xLim;
    case 'tight'
        xLim = xlim; if xLim(1)>0; xLim(1) = 0; end; if xLim(2)<0; xLim(2) = 0; end;
        yLim = ylim; if yLim(1)>0; yLim(1) = 0; end; if yLim(2)<0; yLim(2) = 0; end;
    otherwise
        dbstack; error('Invalid window type')
end
dLim = max(diff(xLim),diff(yLim)).*0.03;
xLim = xLim + [-dLim dLim];
yLim = yLim + [-dLim dLim];
set(ax,'XLim',xLim,'YLim',yLim);
hL(1) = line(xLim,[0 0],'Color','w');
hL(2) = line([0 0],yLim,'Color','w');
uistack(hL,'bottom');
xlabel('real')
ylabel('imag')
ax.Color = 'k';

switch windowType
    case 'full'
    case 'tight'
        grid on
        dTick = min(mean(diff(ax.XTick)),mean(diff(ax.YTick)));
        xTick = round(xLim./dTick).*dTick; xTick = linspace(xTick(1),xTick(2),range(xTick)/dTick+1);
        yTick = round(yLim./dTick).*dTick; yTick = linspace(yTick(1),yTick(2),range(yTick)/dTick+1);
        set(ax,'XTick',xTick,'YTick',yTick);
    otherwise
        dbstack; error('Invalid window type')
end

theta = linspace(0,2*pi,100);
M = sqrt(sum(mean(d{dVencLabel==inf},1).^2));
x = M*cos(theta);
y = M*sin(theta);
hC = plot(ax,x,y,'w');
uistack(hC,'bottom');

% cb = colorbar;
% cb.Colormap = cMap;
% cb.Limits = [0 1];
% cb.Ticks = cMapM1;
% cb.TickLabels = strcat('\pi/(',replace(cellstr(num2str(dVencLabel)),' ',''),'\gamma*100)');
% ylabel(cb,'velocity encoding strenght (T*s^2/m)');
% % get(cb)


% title(legend(hP,replace(cellstr(num2str(dVencLabel)),' ',''),'Location','bestOutside','AutoUpdate','off'),[num2str((1-prob(1).alpha)*100) '%CI of mean' newline 'V_{enc} (cm/s)']);
title(legend(hP,replace(cellstr(num2str(dVencLabel)),' ',''),'Location','bestOutside','AutoUpdate','off'),['V_{enc} (cm/s)']);
uistack(hP,'top');

hold off




% if exist('fOut','var') && ~isempty(fOut)
%     drawnow;
%     saveas(f,[fOut '.png']);
%     saveas(f,[fOut '.fig']);
%     saveas(f,[fOut '.eps']);
% end





