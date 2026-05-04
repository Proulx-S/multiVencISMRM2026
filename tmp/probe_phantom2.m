addpath(genpath('/scratch/bass/tools/util'));
addpath('/scratch/bass/projects/multiVencISMRM2026');
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
% DO NOT conjugate yet

fprintf('size(data): '); disp(size(data));
fprintf('isreal(data): %d\n', isreal(data));
fprintf('dataVenc has Inf: %d\n', any(dataVenc==inf));

% Raw phase at center pixel (12,9) for bestVenc vs inf
[nFE, nPE, nMeas] = size(data);
centerFE = round(nFE/2); centerPE = round(nPE/2);

infIdx  = find(dataVenc==inf);
bestIdx = find(dataVenc==9);
fprintf('n(venc==inf): %d, n(venc==9): %d\n', length(infIdx), length(bestIdx));

if ~isempty(bestIdx)
    px_inf  = squeeze(data(centerFE, centerPE, infIdx(1)));
    px_best = squeeze(data(centerFE, centerPE, bestIdx(1)));
    fprintf('pixel at (center,inf) = %.4g+%.4gi  |abs|=%.4g  phase=%.4f rad\n', real(px_inf), imag(px_inf), abs(px_inf), angle(px_inf));
    fprintf('pixel at (center,v=9) = %.4g+%.4gi  |abs|=%.4g  phase=%.4f rad\n', real(px_best), imag(px_best), abs(px_best), angle(px_best));
    PD_raw = angle(px_best) - angle(px_inf);
    fprintf('Phase difference (raw) = %.4f rad  (vel = %.3f cm/s)\n', PD_raw, PD_raw/pi*9);
end
