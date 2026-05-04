addpath(genpath('/scratch/bass/tools/util'));
phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
data = conj(data);

ID = 6.35; OD = 11.11; bestVenc = 9;

% Check data at center pixel (14,9) manually
cFE = 14; cPE = 9;
bestIdx = squeeze(dataVenc==bestVenc);
infIdx  = squeeze(dataVenc==inf);

fprintf('Number of bestVenc scans: %d, Number of inf scans: %d\n', sum(bestIdx), sum(infIdx));

pixel_inf  = data(cFE, cPE, infIdx);
pixel_best = data(cFE, cPE, bestIdx);
fprintf('pixel_inf phases (after conj): ');
disp(angle(pixel_inf(:)'));
fprintf('pixel_best phases (after conj): ');
disp(angle(pixel_best(:)'));

mean_inf = mean(data(cFE,cPE,infIdx), 3);
fprintf('mean(inf): %.4g+%.4gi, phase=%.4f\n', real(mean_inf), imag(mean_inf), angle(mean_inf));

ref_phase_exp = exp(1j.*angle(mean_inf));
fprintf('ref_phase_exp: %.4g+%.4gi\n', real(ref_phase_exp), imag(ref_phase_exp));

best_normalized = data(cFE,cPE,bestIdx) ./ ref_phase_exp;
fprintf('best_normalized phases: ');
disp(angle(squeeze(best_normalized(:))));
PD_center = angle(mean(best_normalized, 3));
fprintf('PD at center: %.4f rad, vel=%.4f cm/s\n', PD_center, PD_center/pi*bestVenc);

% Now compute full PD map
ref_map = exp(1j.*angle(mean(data(:,:,infIdx),3)));
fprintf('size(ref_map): '); disp(size(ref_map));
fprintf('size(data(:,:,bestIdx)): '); disp(size(data(:,:,bestIdx)));
best_all = data(:,:,bestIdx) ./ ref_map;
fprintf('size(best_all): '); disp(size(best_all));
PD_map = angle(mean(best_all, 3));
fprintf('size(PD_map): '); disp(size(PD_map));
fprintf('PD_map(14,9)=%.4f rad, vel=%.4f cm/s\n', PD_map(cFE,cPE), PD_map(cFE,cPE)/pi*bestVenc);
fprintf('PD_map range: [%.4f, %.4f] rad\n', min(PD_map(:)), max(PD_map(:)));
