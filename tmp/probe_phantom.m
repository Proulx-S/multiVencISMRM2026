addpath(genpath('/scratch/bass/tools/util'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

phantom03dataFile = '/scratch/bass/projects/multiVencISMRM2026/tmp/phantom03.mat';
load(phantom03dataFile, 'data', 'dataVenc', 'dataRun', 'dataNoFlow', 'PEspacing', 'FEspacing');
data = conj(data); dataNoFlow = conj(dataNoFlow);

ID = 6.35; OD = 11.11; bestVenc = 9;

M0 = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos = linspace(FEspacing/2, size(M0,1)*FEspacing-FEspacing/2, size(M0,1));
PEpos = linspace(PEspacing/2, size(M0,2)*PEspacing-PEspacing/2, size(M0,2));
[FEgrid, PEgrid] = ndgrid(FEpos, PEpos);
total = sum(M0(:));
com(1) = sum(FEgrid(:).*M0(:))/total;
com(2) = sum(PEgrid(:).*M0(:))/total;
FEgrid = FEgrid-com(1); PEgrid = PEgrid-com(2);
rGrid = sqrt(PEgrid.^2+FEgrid.^2);
dFE = abs(FEgrid); dPE = abs(PEgrid);
d_far  = sqrt((dFE+FEspacing/2).^2+(dPE+PEspacing/2).^2);
d_near = sqrt(max(0,dFE-FEspacing/2).^2+max(0,dPE-PEspacing/2).^2);
maskBloodOnly = d_far < ID/2;

fprintf('maskBloodOnly: %d pixels\n', sum(maskBloodOnly(:)));
fprintf('size(data): '); disp(size(data));
fprintf('unique(dataVenc): '); disp(unique(dataVenc)');
fprintf('bestVenc in dataVenc: %d\n', any(dataVenc==bestVenc));

PD = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))), 3));
velPD = phase2vel(PD, bestVenc);

fprintf('mean(velPD(maskBloodOnly)) = %.4f cm/s\n', mean(velPD(maskBloodOnly)));
fprintf('range velPD blood: [%.4f, %.4f] cm/s\n', min(velPD(maskBloodOnly)), max(velPD(maskBloodOnly)));
fprintf('mean(rGrid(maskBloodOnly)) = %.4f mm\n', mean(rGrid(maskBloodOnly)));
