addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));

[voxG, spG, nPerVox] = setGrid(4*6.35, 5*4.0, 4, 5, (2^8+1)^2);

fprintf('voxGrid: fovFE=%.2f fovPE=%.2f  matFE=%d matPE=%d  dFE=%.4f dPE=%.4f mm\n', ...
    voxG.fovFE, voxG.fovPE, voxG.matFE, voxG.matPE, voxG.dFE, voxG.dPE);
fprintf('voxGrid.coorFE [mm]: '); fprintf('%.4f  ', voxG.coorFE); fprintf('\n');
fprintf('voxGrid.coorPE [mm]: '); fprintf('%.4f  ', voxG.coorPE); fprintf('\n');
fprintf('voxel boundaries FE: '); fprintf('%.4f  ', voxG.coorFE(1:end-1)+voxG.dFE/2); fprintf('mm\n');
fprintf('voxel boundaries PE: '); fprintf('%.4f  ', voxG.coorPE(1:end-1)+voxG.dPE/2); fprintf('mm\n');
fprintf('\n');
fprintf('spinGrid: matFE=%d matPE=%d  dFE=%.4f dPE=%.4f mm  nPerVox=%d  nTotal=%d\n', ...
    spG.matFE, spG.matPE, spG.dFE, spG.dPE, nPerVox, spG.matFE*spG.matPE);
fprintf('\n');

voxIdx = getVoxIdx(voxG, spG);
for idx = 0:5
    mask = voxIdx == idx;
    cFE = mean(spG.coorFE(any(mask,1)));
    cPE = mean(spG.coorPE(any(mask,2)));
    fprintf('voxel %d: center=(FE=%.3f, PE=%.3f) mm  dist=%.3f mm  nSpins=%d\n', ...
        idx, cFE, cPE, sqrt(cFE^2+cPE^2), sum(mask(:)));
end
