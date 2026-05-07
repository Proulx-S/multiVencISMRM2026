% testWallMask.m
% Experiment with gray overlays for wall voxels on phantom velocity map.
% The maskWallLowMag mask identifies wall voxels (low magnitude pixels
% between the inner and outer vessel walls).

PROJ = '/scratch/bass/projects/multiVencISMRM2026';
addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));

% ── Load data ────────────────────────────────────────────────────────
load(fullfile(PROJ,'tmp','phantom03.mat'), 'data','dataVenc','PEspacing','FEspacing');
data = conj(data);
ID = 6.35; OD = 11.11; bestVenc = 10;

% Coordinates centred on vessel
M      = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos  = linspace(FEspacing/2, size(M,1)*FEspacing-FEspacing/2, size(M,1));
PEpos  = linspace(PEspacing/2, size(M,2)*PEspacing-PEspacing/2, size(M,2));
[FEgrid,PEgrid] = ndgrid(FEpos,PEpos);
com    = [sum(FEgrid(:).*M(:)) sum(PEgrid(:).*M(:))] / sum(M(:));
FEpos  = FEpos - com(1);  PEpos = PEpos - com(2);
[FEgrid,PEgrid] = ndgrid(FEpos,PEpos);

% Masks and velocity map
[~,~,~,~,maskWallLowMag] = makeVesselMasks(FEgrid,PEgrid,FEspacing,PEspacing,ID,OD,M);
velMap = phase2vel(angle(mean(data(:,:,dataVenc==bestVenc),3)), vencToM1(bestVenc));
vLim   = max(abs(velMap(:)));

% Gray RGB overlay layer (neutral mid-gray)
grayRGB  = 0.45 * ones(size(velMap,1), size(velMap,2), 3);
mask     = double(maskWallLowMag);

% ── Figure: 4 overlay styles ────────────────────────────────────────
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 32 16]);
hT = tiledlayout(f,1,4,'TileSpacing','compact','Padding','compact');

styles = {
    'Opaque gray',        1.0;
    'Semi-transparent',   0.6;
    'Light veil',         0.3;
    'Stipple (alpha~0.5)',0.5;
};

for k = 1:4
    ax = nexttile(hT);
    imagesc(ax, PEpos, FEpos, velMap, [-vLim vLim]); axis image;
    ax.Colormap = blueBlackRed; set(ax,'XTick',[],'YTick',[]);
    hold(ax,'on');

    switch k
        case {1,2,3}
            % Solid gray overlay with varying alpha
            alpha = styles{k,2};
            h = image(ax, PEpos, FEpos, grayRGB);
            h.AlphaData = mask * alpha;

        case 4
            % Stippled: checkerboard alpha pattern within mask
            [rows,cols] = size(velMap);
            checker = mod(meshgrid(1:cols,1:rows) + meshgrid(1:rows,1:cols)', 2);
            h = image(ax, PEpos, FEpos, grayRGB);
            h.AlphaData = mask .* (0.3 + 0.5*checker);
    end

    title(ax, styles{k,1});
    colorbar(ax);
end

% ── Harmonise CLims and save ─────────────────────────────────────────
outPath = fullfile(PROJ,'figures','testWallMask.png');
set(findall(f,'Type','axes'),'FontSize',12);
exportgraphics(f, outPath, 'Resolution',150);
fprintf('Saved: %s\n', outPath);
