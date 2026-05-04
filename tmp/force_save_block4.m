% Force regeneration of phantomSummary with simulation overlay
addpath(genpath('/scratch/bass/tools/util'));
addpath(genpath('/scratch/bass/tools/pcMRAsim'));
addpath(genpath('/scratch/bass/tools/red-blue-colormap'));
addpath('/scratch/bass/projects/multiVencISMRM2026');

info.project.storage = '/local/users/Proulx-S/projects/multiVencISMRM2026';
info.project.scratch = '/scratch/bass/projects/multiVencISMRM2026/tmp';
info.project.code    = '/scratch/bass/projects/multiVencISMRM2026';

load(fullfile(info.project.scratch,'phantom03.mat'),'data','dataVenc','dataRun','dataNoFlow','PEspacing','FEspacing');
data = conj(data); dataNoFlow = conj(dataNoFlow);
ID=6.35; OD=11.11; bestVenc=9;

M0=squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
FEpos=linspace(FEspacing/2,size(M0,1)*FEspacing-FEspacing/2,size(M0,1));
PEpos=linspace(PEspacing/2,size(M0,2)*PEspacing-PEspacing/2,size(M0,2));
[FEgrid,PEgrid]=ndgrid(FEpos,PEpos);
total=sum(M0(:)); com(1)=sum(FEgrid(:).*M0(:))/total; com(2)=sum(PEgrid(:).*M0(:))/total;
FEgrid=FEgrid-com(1); FEpos=FEpos-com(1); PEgrid=PEgrid-com(2); PEpos=PEpos-com(2);
rGrid=sqrt(PEgrid.^2+FEgrid.^2);
dFE=abs(FEgrid); dPE=abs(PEgrid);
d_far=sqrt((dFE+FEspacing/2).^2+(dPE+PEspacing/2).^2);
d_near=sqrt(max(0,dFE-FEspacing/2).^2+max(0,dPE-PEspacing/2).^2);
maskBloodOnly=d_far<ID/2; maskTissueOnly=d_near>OD/2;
maskWallLowMag=single(M0<0.44e-7);
theta=linspace(0,2*pi,360);

% Radial profile fit
PD_tmp=angle(mean(data(:,:,dataVenc==bestVenc)./exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))),3));
m1best=vencToM1(bestVenc); velPD_tmp=phase2vel(PD_tmp,m1best);
if mean(velPD_tmp(maskBloodOnly))<0; velPD_tmp=-velPD_tmp; end
rFit3=double(rGrid(maskBloodOnly)); vFit3=double(velPD_tmp(maskBloodOnly));
parabola_fun3=@(vMax,r) vMax.*(1-(r./(ID/2)).^2);
opts3=optimoptions('lsqcurvefit','Display','off');
vMax_fit3=lsqcurvefit(parabola_fun3,max(vFit3),rFit3(:),vFit3(:),0,[],opts3);
vMean_fit3=vMax_fit3/2;
fprintf('vMax=%.3f vMean=%.3f\n', vMax_fit3, vMean_fit3);

saveThis=1;

% Re-run the full block 4 from doIt.m (copy here for force-save)
f = figure('MenuBar','none','ToolBar','none','Units','centimeters','Position',[0 0 35 18.5]);
hT = tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor'); ax = {};

ax{end+1} = nexttile;
M = squeeze(abs(mean(data(:,:,dataVenc==inf),3)));
hIm = imagesc(ax{end},PEpos,FEpos,M,[0 max(M(:))]); axis image;
ylabel(colorbar('Location','westoutside'), 'MR magn. [a.u.]');
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[]); title(ax{end},'phantom ROI');

ax{end+1} = nexttile;
PD   = angle(mean(data(:,:,dataVenc==bestVenc) ./ exp(1j.*angle(mean(data(:,:,dataVenc==inf),3))),3));
CD   = mean(data(:,:,dataVenc==bestVenc),3)-mean(data(:,:,dataVenc==inf),3);
[velCD,phi,velPD,~,~,~] = getPlugFlowEstimates(bestVenc,CD,[],[],PD);
hIm = imagesc(ax{end},PEpos,FEpos,velPD,[-max(abs(velPD(:))) max(abs(velPD(:)))]); axis image;
ylabel(colorbar('Location','westoutside'), 'velocity [cm/s]');
ax{end}.Colormap = redblue; set(ax{end},'XTick',[],'YTick',[]);
title(ax{end},['venc=' num2str(bestVenc) 'cm/s']);

ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskBloodOnly,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[],'Color','none'); hold(ax{end},'on');
plot(ax{end},ID/2*cos(theta),ID/2*sin(theta),'m'); title(ax{end},'blood-only mask');
ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskTissueOnly,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[],'Color','none'); hold(ax{end},'on');
plot(ax{end},OD/2*cos(theta),OD/2*sin(theta),'m'); title(ax{end},'tissue-only mask');
ax{end+1} = nexttile;
imagesc(ax{end},PEpos,FEpos,maskWallLowMag,[0 1]); axis image;
ax{end}.Colormap = gray; set(ax{end},'XTick',[],'YTick',[],'Color','none'); hold(ax{end},'on');
plot(ax{end},ID/2*cos(theta),ID/2*sin(theta),'m');
plot(ax{end},OD/2*cos(theta),OD/2*sin(theta),'m'); title(ax{end},'low-mag wall mask');

% Run matched simulation
pDefault4 = runSim;
pVessel4 = pDefault4.pVessel; pVessel4.ID=ID; pVessel4.WT=2.38125;
pVessel4.vMean=vMean_fit3; pVessel4.profile='parabolic1'; pVessel4.S.lumen=[];
pSim4 = pDefault4.pSim; pSim4.fovFE=FEspacing; pSim4.fovPE=PEspacing;
pSim4.matFE=3; pSim4.matPE=3; pSim4.nSpin=(2^8+1)^2;
pMri4 = pDefault4.pMri; pMri4.venc.method='PCmono';
pMri4.venc.vencList=sort(unique(dataVenc(~isinf(dataVenc))),'descend');
pMri4.venc.FVEres=0; pMri4.venc.FVEbw=0; pMri4.venc.FVEvel=[]; pMri4.venc.vencMin=[]; pMri4.venc.vencMax=[];
res4 = runSim(pVessel4, pSim4, pMri4);
vencListSim4 = sort(unique(pMri4.venc.vencList),'descend');
Iref4=res4.I(1,1,1,1,1,2); Ienc4=squeeze(res4.I(1,1,1,1,:,1));
Isim4_all=[Iref4;Ienc4(:)]; vsim4_all=[inf;vencListSim4(:)];

ax{end+1} = nexttile([3 3]);
I=squeeze(mean(data,[1 2])); Ivenc=squeeze(dataVenc); Irun=squeeze(dataRun);
plotMultiVenc(ax{end},I,Ivenc,Irun,[],'hot');
hold(ax{end},'on');
Mnorm4=abs(Iref4);
cMapSim4=cool(length(vencListSim4));
for kk4=1:length(vencListSim4)
    v4_=vencListSim4(kk4);
    plot(ax{end}, real(Isim4_all(vsim4_all==v4_))/Mnorm4, imag(Isim4_all(vsim4_all==v4_))/Mnorm4, ...
         's','MarkerFaceColor',cMapSim4(kk4,:),'MarkerEdgeColor','w','MarkerSize',5);
end

drawnow;
saveas(f, fullfile(info.project.storage,'figures','phantomSummary.fig'));
exportgraphics(f, fullfile(info.project.storage,'figures','phantomSummary.png'));
exportgraphics(f, fullfile(info.project.storage,'figures','phantomSummary.svg'));
disp('Saved phantomSummary with simulation overlay');
