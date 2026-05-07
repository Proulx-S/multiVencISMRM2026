% testColormap.m — test perceptual red/blue diverging colormap on phantom data

addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));

% Ensure colorspace conversion tool is available (same dependency as colormap_bivariateBlackToSpectral)
if exist('colorspace','file') ~= 2
    toolDir = '/scratch/bass/tools';
    tool    = 'Colorspace-Transformations';
    toolURL = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/28790/versions/5/download/zip';
    if ~exist(fullfile(toolDir,tool),'dir')
        tmpZip = fullfile(tempdir,'colorspace.zip');
        websave(tmpZip, toolURL); unzip(tmpZip, fullfile(toolDir,tool)); delete(tmpZip);
    end
    addpath(genpath(fullfile(toolDir,tool)));
end

% ── Load phantom data ─────────────────────────────────────────────────
load('/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat', ...
     'data','dataVenc','PEspacing','FEspacing');
data = conj(data);

% Centre coordinates on vessel centre of mass
M      = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos  = linspace(FEspacing/2, size(M,1)*FEspacing-FEspacing/2, size(M,1));
PEpos  = linspace(PEspacing/2, size(M,2)*PEspacing-PEspacing/2, size(M,2));
[FEg,PEg] = ndgrid(FEpos,PEpos);
com    = [sum(FEg(:).*M(:)) sum(PEg(:).*M(:))] / sum(M(:));
FEpos  = FEpos - com(1);
PEpos  = PEpos - com(2);

% Phase map at bestVenc (good test: vessel + background, symmetric range)
bestVenc = 10;
phaseMap = angle(mean(data(:,:,dataVenc==bestVenc),3));


% ── Build perceptual colormap in CIE L*C*H° space ────────────────────
% Control points: [t, L, C, H_deg]
%   CIE LCH hue convention: 0°=red, 90°=yellow, 180°=green, 270°=blue
%
%   light blue  →  pure blue  →  darker blue  →  gray  →  dark red  →  pure red  →  light red
ctrl = [
    0.00,  80,  25,  280;   % light blue   (sky blue:  high L, low-mid C)
    0.28,  38,  50,  280;   % pure blue    (saturated: low-mid L, high C)
    0.40,  12,  18,  285;   % darker blue  (deep navy: low L, mid C)
    0.50,   0,   0,    0;   % black        (achromatic: L=0, C=0)
    0.60,  12,  18,   30;   % dark red     (deep crimson: low L, mid C)
    0.72,  42,  60,   25;   % pure red     (saturated: low-mid L, high C)
    1.00,  80,  28,   18;   % light red    (pink-red: high L, low-mid C)
];

N   = 256;
t   = linspace(0, 1, N)';
L   = interp1(ctrl(:,1), ctrl(:,2), t, 'pchip');
C   = interp1(ctrl(:,1), ctrl(:,3), t, 'pchip');
H   = interp1(ctrl(:,1), ctrl(:,4), t, 'pchip');
C   = max(C, 0);    % chroma must be non-negative

% Convert to RGB via CIE LCH → sRGB
LCH_img  = reshape([L, C, H], [N, 1, 3]);         % N×1×3
LCH_img  = permute(LCH_img, [2 1 3]);             % 1×N×3
RGB_raw  = colorspace('LCH->RGB', LCH_img);       % 1×N×3
cmap     = squeeze(max(0, min(1, RGB_raw)));       % N×3, clipped to sRGB gamut


% ── Figure: comparison ───────────────────────────────────────────────
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 36 14]);
hT = tiledlayout(f, 2, 3, 'TileSpacing','compact','Padding','compact');

% Row 1: phase maps
ax = nexttile(hT, 1);
imagesc(ax, PEpos, FEpos, phaseMap, [-pi pi]); axis image;
ax.Colormap = redblue; set(ax,'XTick',[],'YTick',[]);
ylabel(colorbar(ax), 'phase [rad]'); title(ax, 'redblue (reference)');

ax = nexttile(hT, 2);
imagesc(ax, PEpos, FEpos, phaseMap, [-pi pi]); axis image;
ax.Colormap = cmap; set(ax,'XTick',[],'YTick',[]);
ylabel(colorbar(ax), 'phase [rad]'); title(ax, 'perceptual (new)');

% Row 1 col 3: colormap strips side by side
ax = nexttile(hT, 3);
strip = cat(1, permute(redblue(N),[3 1 2]), permute(cmap,[3 1 2]));  % 2×N×3
strip = permute(strip,[1 2 3]);
image(ax, linspace(-pi,pi,N), [1 2], strip);
set(ax,'YTick',[1 2],'YTickLabel',{'redblue','perceptual'},'XAxisLocation','bottom');
xlabel(ax,'phase [rad]'); title(ax,'colormap strips');

% Row 2: luminance profiles (perceptual smoothness check)
ax = nexttile(hT, 4);
cmapLAB = colorspace('RGB->Lab', permute(redblue(N),[3 1 2]));
plot(ax, linspace(-pi,pi,N), squeeze(cmapLAB(1,:,1)), 'k', 'LineWidth',1.5);
xlabel(ax,'phase [rad]'); ylabel(ax,'L* (perceptual lightness)');
title(ax,'redblue — L* profile'); ylim(ax,[0 100]); grid(ax,'on');

ax = nexttile(hT, 5);
cmapLAB = colorspace('RGB->Lab', permute(cmap,[3 1 2]));
plot(ax, linspace(-pi,pi,N), squeeze(cmapLAB(1,:,1)), 'k', 'LineWidth',1.5);
xlabel(ax,'phase [rad]'); ylabel(ax,'L* (perceptual lightness)');
title(ax,'perceptual — L* profile'); ylim(ax,[0 100]); grid(ax,'on');

% Row 2 col 3: chroma profile
ax = nexttile(hT, 6);
cmapLAB_rb  = colorspace('RGB->Lab', permute(redblue(N),[3 1 2]));
cmapLAB_new = colorspace('RGB->Lab', permute(cmap,     [3 1 2]));
C_rb  = squeeze(sqrt(cmapLAB_rb( 1,:,2).^2 + cmapLAB_rb( 1,:,3).^2));
C_new = squeeze(sqrt(cmapLAB_new(1,:,2).^2 + cmapLAB_new(1,:,3).^2));
ph = linspace(-pi,pi,N);
plot(ax, ph, C_rb,  'b', 'LineWidth',1.5); hold(ax,'on');
plot(ax, ph, C_new, 'r', 'LineWidth',1.5);
legend(ax,{'redblue','perceptual'},'Location','north');
xlabel(ax,'phase [rad]'); ylabel(ax,'C* (chroma)');
title(ax,'chroma profiles'); ylim(ax,[0 100]); grid(ax,'on');

% ── Save ─────────────────────────────────────────────────────────────
outDir  = '/scratch/bass/projects/multiVencISMRM2026/figures';
if ~exist(outDir,'dir'); mkdir(outDir); end
outPath = fullfile(outDir, 'testColormap.png');
exportgraphics(f, outPath, 'Resolution', 150);
fprintf('Saved: %s\n', outPath);
