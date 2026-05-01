function [data, dataVenc,PEspacing,FEspacing] = loadPhantom03(dataPath)
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

imgAllFlow   = [];
imgAllNoFlow = [];
PDflow       = [];
PDnoFlow     = [];
vencAllNoFlow = [];
vencAllFlow   = [];
vencPDflow    = [];
vencPDnoFlow  = [];
runAllNoFlow = [];
runAllFlow   = [];
runPDflow    = [];
runPDnoFlow  = [];
for iRun = 1:length(runList)
    load(fullfile(runList(iRun).folder,runList(iRun).name));
    img = conj(img); % flip sign of phase
    if size(venc,11)~=size(img,11) && size(venc,11)==1
        venc = repmat(venc,[1 1 1 1 1 1 1 1 1 1 size(img,11) 1 1 1 1 1]);
    end
    run = ones(size(venc)).*iRun;

    % Remove coil phase
    if size(img,4)>1
        img = mean(img .* exp(-1i.*angle(mean(kCoil,11))),4);
    end
    
    % Compute PD
    refInd = squeeze(              venc(:,:,:,:,:,:,     :,:,:,:,1,:,:,:,:,:)==inf);
    PD     = img .* exp(-1i * angle(img(:,:,:,:,:,:,refInd,:,:,:,1,:,:,:,:,:))    );

    
    % Compile for signal averaging -- time points with no flow
    stillIdx2  = ((  timeOffset(iRun) + 2 + infuse + withdraw  ):cycleLength:(size(img,11)+cycleLength));
    stillIdx{iRun}  = sort([stillIdx2(:)])';
    if stillIdx{iRun}(1)>cycleLength
        stillIdx{iRun} = stillIdx{iRun} - cycleLength;
    end
    stillIdx{iRun}(stillIdx{iRun}>size(img,11)) = [];

    infuseIdx2 = ((  timeOffset(iRun) + 2                      ):cycleLength:(size(img,11)+cycleLength));
    infuseIdx{iRun} = sort([infuseIdx2(:)])';
    if infuseIdx{iRun}(1)>cycleLength
        infuseIdx{iRun} = infuseIdx{iRun} - cycleLength;
    end
    infuseIdx{iRun}(infuseIdx{iRun}>size(img,11)) = [];

    tmp = img(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    imgAllNoFlow = cat(7,imgAllNoFlow,mean(tmp(:,:,:,:,:,:,:),4));
    tmp = venc(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    vencAllNoFlow      = cat(7,vencAllNoFlow,tmp(:,:,:,:,:,:,:));
    tmp = run(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    runAllNoFlow      = cat(7,runAllNoFlow,tmp(:,:,:,:,:,:,:));
    tmp = PD(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    PDnoFlow = cat(7,PDnoFlow,mean(tmp(:,:,:,:,:,:,:),4));
    tmp = venc(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    vencPDnoFlow = cat(7,vencPDnoFlow,tmp(:,:,:,:,:,:,:));
    tmp = run(:,:,:,:,:,:,:,:,:,:,stillIdx{iRun},:,:,:,:,:);
    runPDnoFlow = cat(7,runPDnoFlow,tmp(:,:,:,:,:,:,:));

    tmp = img(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    imgAllFlow = cat(7,imgAllFlow,mean(tmp(:,:,:,:,:,:,:),4));
    tmp = venc(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    vencAllFlow = cat(7,vencAllFlow,tmp(:,:,:,:,:,:,:));
    tmp = run(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    runAllFlow = cat(7,runAllFlow,tmp(:,:,:,:,:,:,:));
    tmp = PD(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    PDflow = cat(7,PDflow,mean(tmp(:,:,:,:,:,:,:),4));
    tmp = venc(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    vencPDflow = cat(7,vencPDflow,tmp(:,:,:,:,:,:,:));
    tmp = run(:,:,:,:,:,:,:,:,:,:,infuseIdx{iRun},:,:,:,:,:);
    runPDflow = cat(7,runPDflow,tmp(:,:,:,:,:,:,:));
end

% sort vencs
[a,b] = sort(vencAllNoFlow);
imgAllNoFlow = imgAllNoFlow(:,:,:,:,:,:,b);
vencAllNoFlow = vencAllNoFlow(:,:,:,:,:,:,b);
runAllNoFlow = runAllNoFlow(:,:,:,:,:,:,b);

[a,b] = sort(vencAllFlow);
imgAllFlow = imgAllFlow(:,:,:,:,:,:,b);
vencAllFlow = vencAllFlow(:,:,:,:,:,:,b);
runAllFlow = runAllFlow(:,:,:,:,:,:,b);

[a,b] = sort(vencPDflow);
PDflow = PDflow(:,:,:,:,:,:,b);
vencPDflow = vencPDflow(:,:,:,:,:,:,b);
runPDflow = runPDflow(:,:,:,:,:,:,b);

[a,b] = sort(vencPDnoFlow);
PDnoFlow = PDnoFlow(:,:,:,:,:,:,b);
vencPDnoFlow = vencPDnoFlow(:,:,:,:,:,:,b);
runPDnoFlow = runPDnoFlow(:,:,:,:,:,:,b);



vencPDflowList   = sort(unique(vencPDflow));
vencPDnoFlowList = sort(unique(vencPDnoFlow));
runPDflowList   = sort(unique(runPDflow));
runPDnoFlowList = sort(unique(runPDnoFlow));
% average time within each run
dataFlow_tAv   = nan(size(PDflow  ,1),size(PDflow  ,2),length(vencPDflowList  ),length(runPDflowList));
dataNoFlow_tAv = nan(size(PDnoFlow,1),size(PDnoFlow,2),length(vencPDnoFlowList),length(runPDnoFlowList));
for iRun = 1:length(runPDflowList)
    for iVenc = 1:length(vencPDflowList)
        if nnz(vencPDflow==vencPDflowList(iVenc)   &   runPDflow==runPDflowList(iRun))
            dataFlow_tAv(:,:,iVenc,iRun)   = mean(PDflow  (:,:,:,:,:,:,  vencPDflow==vencPDflowList(iVenc)   &   runPDflow==runPDflowList(iRun)  ),7);
            dataNoFlow_tAv(:,:,iVenc,iRun) = mean(PDnoFlow(:,:,:,:,:,:,vencPDnoFlow==vencPDnoFlowList(iVenc) & runPDnoFlow==runPDnoFlowList(iRun)),7);
        end
    end
end
% correct eddy current effects
dataFlow_tAv = dataFlow_tAv .* exp(-1i*angle(dataNoFlow_tAv));
% average over runs
data = mean(permute(dataFlow_tAv,[3 1 2 4]),4,"omitnan"); % [venc x FE x PE], use nanmean to ignore NaNs
dataVenc = vencPDnoFlowList;
    