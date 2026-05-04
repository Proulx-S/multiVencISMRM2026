addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

p = runSim; % defaults
p.pSim.fovFE = 15; p.pSim.fovPE = 18; % mm -- big enough to show full OD=11.11mm vessel
p.pVessel.ID = 6.35; p.pVessel.WT = 2.38125; p.pVessel.vMean = 4.234;
p.pMri.venc.method = 'PCmono';
p.pMri.venc.vencList = [100; 10]; % two venc, just to test
p.pMri.venc.FVEres=0; p.pMri.venc.FVEbw=0; p.pMri.venc.FVEvel=[]; p.pMri.venc.vencMin=[]; p.pMri.venc.vencMax=[];
res = runSim(p.pVessel, p.pSim, p.pMri, [], false);

fprintf('Fields in res: '); disp(fieldnames(res)');
fprintf('Fields in res.pSim: '); disp(fieldnames(res.pSim)');
fprintf('Fields in res.pVessel: '); disp(fieldnames(res.pVessel)');
fprintf('Fields in res.pVessel.mask: '); disp(fieldnames(res.pVessel.mask)');
fprintf('size(res.magMap)='); disp(size(res.magMap));
fprintf('size(res.vMap)='); disp(size(res.vMap));
fprintf('size(res.pVessel.mask.lumen)='); disp(size(res.pVessel.mask.lumen));
fprintf('pSim fields with grid: '); disp(fieldnames(res.pSim)');
if isfield(res.pSim,'gridFE'); fprintf('gridFE size: '); disp(size(res.pSim.gridFE)); end
if isfield(res.pSim,'gridPE'); fprintf('gridPE size: '); disp(size(res.pSim.gridPE)); end
