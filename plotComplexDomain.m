function ax = plotComplexDomain(ax, I, windowType, plotStyle)
% plotComplexDomain  Complex-plane signal plot with unit circle and crosshairs.
%
%   ax = plotComplexDomain(ax, I, windowType, plotStyle)
%
%   ax         : target axes handle ([] → use current axes)
%   I          : complex vector, pre-normalized (caller's responsibility)
%   windowType : 'tight' (default) — origin included, tight to data
%                'full'            — square, origin-centered, equal axes
%   plotStyle  : 'markers' (default) — one dot per entry
%                'line'             — connected line in array order

if isempty(ax);                                        ax = gca;        end
if ~exist('windowType','var') || isempty(windowType);  windowType = 'tight';   end
if ~exist('plotStyle','var')  || isempty(plotStyle);   plotStyle  = 'markers'; end

hold(ax, 'on');

x = real(I(:));
y = imag(I(:));

switch plotStyle
    case 'markers'
        plot(ax, x, y, 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor','w', 'MarkerSize',8, 'LineWidth',1);
    case 'line'
        plot(ax, x, y, 'w-', 'LineWidth', 1.5);
    otherwise
        error('plotStyle must be ''markers'' or ''line''');
end

theta = linspace(0, 2*pi, 360);
plot(ax, cos(theta), sin(theta), 'w--', 'LineWidth', 0.8);
xline(ax, 0, 'w', 'LineWidth', 0.5, 'Alpha', 0.4);
yline(ax, 0, 'w', 'LineWidth', 0.5, 'Alpha', 0.4);

axis(ax, 'image', 'tight');
switch windowType
    case 'full'
        lim  = max(abs([xlim(ax) ylim(ax)]));
        xLim = [-lim lim];
        yLim = [-lim lim];
    case 'tight'
        xLim = xlim(ax); if xLim(1)>0; xLim(1)=0; end; if xLim(2)<0; xLim(2)=0; end;
        yLim = ylim(ax); if yLim(1)>0; yLim(1)=0; end; if yLim(2)<0; yLim(2)=0; end;
        dTick = min(mean(diff(ax.XTick)), mean(diff(ax.YTick)));
        xTick = round(xLim./dTick).*dTick; xTick = linspace(xTick(1),xTick(2),range(xTick)/dTick+1);
        yTick = round(yLim./dTick).*dTick; yTick = linspace(yTick(1),yTick(2),range(yTick)/dTick+1);
        set(ax, 'XTick',xTick, 'YTick',yTick);
    otherwise
        error('windowType must be ''tight'' or ''full''');
end
dLim = max(diff(xLim), diff(yLim)) * 0.03;
set(ax, 'XLim',xLim+[-dLim dLim], 'YLim',yLim+[-dLim dLim]);

ax.Color     = 'k';
ax.GridColor = [0.5 0.5 0.5];
grid(ax, 'on');
xlabel(ax, 'real');
ylabel(ax, 'imag');

hold(ax, 'off');
