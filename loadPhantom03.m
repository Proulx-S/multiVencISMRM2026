function [data, dataVenc, dataRun, dataNoFlow, PEspacing, FEspacing] = loadPhantom03(dataPath)
%%%%%%%%%%%%%%%%%
% Load data cropped to include a bit of static agar around the tube

FEspacing = 0.5;      % mm per pixel in FE direction
PEspacing = 0.8929;   % mm per pixel in PE direction
dataAspectRatio = [FEspacing PEspacing 1];

coilMethod = 'bartMap-';
timeOffset = [3 6 0 2 5];
fileName = ['_fft_coilComb-' strrep(coilMethod, '-', '') '_FEcrop150-176_PEcrop83-100'];
% fileName = '_fft_FEcrop150-176_PEcrop83-100';
runList = dir(fullfile(dataPath,'raw',['*' fileName '.mat'])); runList([runList.isdir]) = [];
runList(contains({runList.name},'_longTR')) = [];


infuse   = 3;
withdraw = 3;
still    = 3;
cycleLength = infuse + withdraw + still;
stillIdx  = cell(size(runList));
infuseIdx = cell(size(runList));

for iRun = 1:length(runList)
    load(fullfile(runList(iRun).folder,runList(iRun).name));
    
    % Get time points with flow
    infuseIdx2 = ((  timeOffset(iRun) + 2                      ):cycleLength:(size(img,11)+cycleLength));
    infuseIdx{iRun} = sort([infuseIdx2(:)])';
    if infuseIdx{iRun}(1)>cycleLength
        infuseIdx{iRun} = infuseIdx{iRun} - cycleLength;
    end
    infuseIdx{iRun}(infuseIdx{iRun}>size(img,11)) = [];
    data{iRun}     = img(  :,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    dataVenc{iRun} = repmat(venc,[1 1 1 1 1 1 1            1 1 1 length(infuseIdx{iRun}) 1 1 1 1]);
    dataRun{iRun}  = repmat(iRun,[1 1 1 1 1 1 size(venc,7) 1 1 1 length(infuseIdx{iRun}) 1 1 1 1]);

    % Get time points with no flow
    stillIdx2  = ((  timeOffset(iRun) + 2 + infuse + withdraw  ):cycleLength:(size(img,11)+cycleLength));
    stillIdx{iRun}  = sort([stillIdx2(:)])';
    if stillIdx{iRun}(1)>cycleLength
        stillIdx{iRun} = stillIdx{iRun} - cycleLength;
    end
    stillIdx{iRun}(stillIdx{iRun}>size(img,11)) = [];
    dataNoFlow{iRun}     = img(  :,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    dataNoFlowVenc{iRun} = repmat(venc,[1 1 1 1 1 1 1            1 1 1 length(stillIdx{iRun}) 1 1 1 1]);
    dataNoFlowRun{iRun}  = repmat(iRun,[1 1 1 1 1 1 size(venc,7) 1 1 1 length(stillIdx{iRun}) 1 1 1 1]);

    
    % ECC using no flow data points
    % data{iRun} = data{iRun}./exp(1j.*angle(mean(dataNoFlow{iRun},11))); % voxel-wise ECC (quite noisy)
    % data{iRun} = data{iRun}./exp(1j.*angle(mean(dataNoFlow{iRun},[1 2 11]))); % ROI-wise ECC (more realistic for a pseudo-voxel ROI)
    data{iRun}       = data{iRun}      ./exp(1j.*angle(mean(dataNoFlow{iRun},[2 11]))); % row-wise ECC (less noisy, leverages the directional nature of the background phase error)
    dataNoFlow{iRun} = dataNoFlow{iRun}./exp(1j.*angle(mean(dataNoFlow{iRun},[2 11]))); % row-wise ECC (less noisy, leverages the directional nature of the background phase error)
end



% Compile across runs and sort by venc and runs
data       = cat(11,data{:}    );
dataVenc   = cat(11,dataVenc{:});
dataRun    = cat(11,dataRun{:} );
dataNoFlow = cat(11,dataNoFlow{:}    );

data       = data(:,:,:);
dataVenc   = dataVenc(:,:,:);
dataRun    = dataRun(:,:,:);
dataNoFlow = dataNoFlow(:,:,:);

% Sort vencs
[~,b] = sort(dataVenc,'descend');
data       = data(:,:,b);
dataVenc   = dataVenc(:,:,b);
dataRun    = dataRun(:,:,b);
dataNoFlow = dataNoFlow(:,:,b);
