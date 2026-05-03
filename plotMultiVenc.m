function plotMultiVenc(ax,I,venc,run)


vencList = sort(unique(venc),'descend');
for vencIdx = 1:size(vencList,1)
    d{vencIdx,1}     = [real(I(venc==vencList(vencIdx))) imag(I(venc==vencList(vencIdx)))];
    dVenc{vencIdx,1} = venc(venc==vencList(vencIdx));
    dVencLabel(vencIdx,1) = unique(dVenc{vencIdx,1});
    dRun{vencIdx,1}  = run( venc==vencList(vencIdx));
end
% dProb = credibleMean2d(d);



cMap = jet;
cMap = interp1(linspace(0,1,size(cMap,1))',cMap,linspace(0,1,length(d))');
% for dIdx = 1:length(d)
%     hP(dIdx) = plot(ax,dProb{dIdx}.CIcontour,'FaceColor',cMap(dIdx,:),'EdgeColor','none','FaceAlpha',1); hold on
% end
for dIdx = 1:length(d)
    hP(dIdx) = plot(ax,mean(d{dIdx,1}(:,1),1),mean(d{dIdx,1}(:,2),1),'o','MarkerFaceColor',cMap(dIdx,:),'MarkerEdgeColor','k'); hold on
end
axis image tight
        


axis image tight
xLim = xlim; if xLim(1)>0; xLim(1) = 0; end; if xLim(2)<0; xLim(2) = 0; end;
yLim = ylim; if yLim(1)>0; yLim(1) = 0; end; if yLim(2)<0; yLim(2) = 0; end;
dLim = max(diff(xLim),diff(yLim)).*0.03;
xLim = xLim + [-dLim dLim];
yLim = yLim + [-dLim dLim];
set(ax,'XLim',xLim,'YLim',yLim);
uistack(xline(0,'w'),'bottom');
uistack(yline(0,'w'),'bottom');
grid on
xlabel('real')
ylabel('imag')
ax.Color = 'k';
dTick = min(mean(diff(ax.XTick)),mean(diff(ax.YTick)));
xTick = round(xLim./dTick).*dTick; xTick = linspace(xTick(1),xTick(2),range(xTick)/dTick+1);
yTick = round(yLim./dTick).*dTick; yTick = linspace(yTick(1),yTick(2),range(yTick)/dTick+1);
set(ax,'XTick',xTick,'YTick',yTick);
theta = linspace(0,2*pi,100);
M = sqrt(sum(mean(d{dVencLabel==inf},1).^2));
x = M*cos(theta);
y = M*sin(theta);
uistack(plot(x,y,'w'),'bottom');

% title(legend(hP,replace(cellstr(num2str(dVencLabel)),' ',''),'Location','bestOutside','AutoUpdate','off'),[num2str((1-prob(1).alpha)*100) '%CI of mean' newline 'V_{enc} (cm/s)']);
title(legend(hP,replace(cellstr(num2str(dVencLabel)),' ',''),'Location','bestOutside','AutoUpdate','off'),['V_{enc} (cm/s)']);
uistack(hP,'top');




% if exist('fOut','var') && ~isempty(fOut)
%     drawnow;
%     saveas(f,[fOut '.png']);
%     saveas(f,[fOut '.fig']);
%     saveas(f,[fOut '.eps']);
% end





