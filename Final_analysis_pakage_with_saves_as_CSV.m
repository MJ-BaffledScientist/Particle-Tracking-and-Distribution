clear all;
close all;
clc;
%foldernames
probenms = {'Untransfected','EGFP-C1','GFP-Tubulin','EMTB-3xGFP','Mn1+2+1+2','Dynein 2219','Kif5c'};
pathname = '/Path/';
foldername = 'Folder Name';
%set nm per pixel
lengthScale = 70e3/270;
%set the s per frame
timeScale = 30e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE THIS FOR DIFFERENT MAXIMUM LENGTHS AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%VELOCITIES
Lover = 500;
Velover = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%STATS
cells_meandisps = {};
cells_meanvels = {};
cells_mediandisps = {};
cells_medianvels = {};
cells_meaninstdisps = {};
cells_meaninstvels = {};
cells_medianinstdisps = {};
cells_medianinstvels = {};
cells_loverdisps = {};
cells_lovervels = {};
cells_lovermeandisps = {};
cells_lovermeanvels = {};
cells_lovermediandisps = {};
cells_lovermedianvels = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%START ANALYSIS
for pbnmdex = 1:length(probenms)
    celld = [];
    cellp = [];
    dispscell=[];
    velcells = [];
    dispsforboxplot = {};
    velsforboxplot = {};
    loverdispsforboxplot = {};
    lovervelsforboxplot = {};
    instdispsforboxplot = {};
    instvelsforboxplot = {};
    namesforboxplot = {};
    boxplotcounter = 1;
    %for stats
    temp_meandisps = [];
    temp_meanvels = [];
    temp_mediandisps = [];
    temp_medianvels = [];
    temp_meaninstdisps = [];
    temp_meaninstvels = [];
    temp_medianinstdisps = [];
    temp_medianinstvels = [];
    temp_loverdisps = [];
    temp_lovervels = [];
    temp_lovermeandisps = [];
    temp_lovermeanvels = [];
    temp_lovermediandisps = [];
    temp_lovermedianvels = [];
    %start experiment loops
    for i = 1:3
        %set full path
        fullpath = [pathname foldername num2str(i)];
        disp(fullpath)
        %read in fluorescence data files
        datafiles = dir([fullpath '/' probenms{pbnmdex} '.txt']);
        s = datafiles.name;
        disp(s)
        data = csvread([fullpath '/' s]);
        cellp = [cellp; data];
        %read in track data files
        trdatafiles = dir([fullpath '/' probenms{pbnmdex} '/data*.mat']);
        trfilenames = {};
        for j = 1:10
            trfilenames{j}= trdatafiles(j).name;
        end
        for j = 1:length(data)
            fstr = [num2str(j) '_t1.mat'];
            strcell = strfind(trfilenames,fstr,'ForceCellOutput',true);
            trindex = find(~cellfun('isempty', strcell));
            w = load([fullpath '/' probenms{pbnmdex} '/' trfilenames{trindex}]);
            [vi,vw] = MeanVel(w.tr,lengthScale,timeScale);
            [di,dw] = MeanDisp(w.tr,lengthScale); 
            %get displacement of each track
            disps = [];
            vels = [];
            for k = 1:length(w.tr)
                disps = [disps; lengthScale*sqrt((w.tr{k}(end,1)-w.tr{k}(1,1)).^2 + (w.tr{k}(end,2)-w.tr{k}(1,2)).^2)];
                %vels = [vels;max((lengthScale/(timeScale))*sqrt((w.tr{k}(2:end,1)-w.tr{k}(1:end-1,1)).^2 + (w.tr{k}(2:end,2)-w.tr{k}(1:end-1,2)).^2))];
                vels = [vels;(lengthScale/(timeScale*length(w.tr{k}(:,1))))*sqrt((w.tr{k}(end,1)-w.tr{k}(1,1)).^2 + (w.tr{k}(end,2)-w.tr{k}(1,2)).^2)];
            end
            %add variables for box plot
            dispsforboxplot{boxplotcounter} = disps(disps>=0);
            velsforboxplot{boxplotcounter} = vels(vels>=0);
            loverdispsforboxplot{boxplotcounter} = numel(disps(disps>Lover))/numel(disps);
            lovervelsforboxplot{boxplotcounter} = numel(vels(disps>Lover))/numel(vels);
            instdispsforboxplot{boxplotcounter} = di(di>=0);
            instvelsforboxplot{boxplotcounter} = vi(vi>=0);
            
            namesforboxplot{boxplotcounter} = ['Exp ' num2str(i) ', Cell ' num2str(boxplotcounter)];
            boxplotcounter = boxplotcounter+1;
            %length over in meters
            dispover = numel(disps(disps>Lover))/numel(disps);
            velover = numel(vels(disps>Lover))/numel(vels);
            %velocity over in meters per second
            velcells = [velcells;velover];
            dispscell = [dispscell;dispover];
            celld = [celld,[median(vi);median(vw);median(di);median(dw)]];
%             disp(velover)
%             disp(numel(vels))
            %for stats
            temp_meandisps = [temp_meandisps; mean(disps)];
            temp_meanvels = [temp_meanvels; mean(vels)];
            temp_mediandisps = [temp_mediandisps; median(disps)];
            temp_medianvels = [temp_medianvels; median(vels)];
            temp_meaninstdisps = [temp_meaninstdisps; mean(di)];
            temp_meaninstvels = [temp_meaninstvels; mean(vi)];
            temp_medianinstdisps = [temp_medianinstdisps; median(di)];
            temp_medianinstvels = [temp_medianinstvels; median(vi)];
            temp_loverdisps = [temp_loverdisps; dispover];
            temp_lovervels = [temp_lovervels; velover];
            temp_lovermeandisps = [temp_lovermeandisps; mean(disps(disps>Lover))];
            temp_lovermeanvels = [temp_lovermeanvels; mean(vels(disps>Lover))];
            temp_lovermediandisps = [temp_lovermediandisps; median(disps(disps>Lover))];
            temp_lovermedianvels = [temp_lovermedianvels; median(vels(disps>Lover))];
        end
    end
    %normalize fluorescence intensities
    fluornorm = cellp(:,1)./cellp(:,2);
    [sorted,ranks]=sort(fluornorm,'ascend');
    %sort cell
    sorteddispsforboxplot={};
    sortedvelsforboxplot={};
    sortedloverdispsforboxplot=zeros(1,length(loverdispsforboxplot));
    sortedlovervelsforboxplot=zeros(1,length(lovervelsforboxplot));
    sortedinstdispsforboxplot={};
    sortedinstvelsforboxplot={};
    sortednamesforboxplot={};
    for rankdex = 1:length(ranks)
        sorteddispsforboxplot{rankdex} = dispsforboxplot{ranks(rankdex)};
        sortedvelsforboxplot{rankdex}=velsforboxplot{ranks(rankdex)};
        sortedloverdispsforboxplot(rankdex)=loverdispsforboxplot{ranks(rankdex)};
        sortedlovervelsforboxplot(rankdex)=lovervelsforboxplot{ranks(rankdex)};
        sortedinstdispsforboxplot{rankdex}=instdispsforboxplot{ranks(rankdex)};
        sortedinstvelsforboxplot{rankdex}=instvelsforboxplot{ranks(rankdex)};
        sortednamesforboxplot{rankdex} = namesforboxplot{ranks(rankdex)};
    end
    %boxplots for each probe
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
    figure
    boxplot2(sorteddispsforboxplot,'labels',sortednamesforboxplot,'orientation','horizontal','whisker',1,'PlotStyle','traditional')
    title([probenms{pbnmdex} ' Displacements Ranked by Fluorescence Intensity (Ascending)'])
    xlabel('Displacements $(nm)$','interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_disps_ranked_boxplot'])
    saveas(gcf,[probenms{pbnmdex} '_disps_ranked_boxplot.png'])
    
    figure
    boxplot2(sortedvelsforboxplot,'labels',sortednamesforboxplot,'orientation','horizontal','whisker',1,'PlotStyle','traditional')
    title([probenms{pbnmdex} ' Velocities Ranked by Fluorescence Intensity (Ascending)'])
    xlabel('Velocities $(nm s^{-1})$','interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_vel_ranked_boxplot'])
    saveas(gcf,[probenms{pbnmdex} '_vel_ranked_boxplot.png'])
    
    figure
    boxplot2(sortedinstdispsforboxplot,'labels',sortednamesforboxplot,'orientation','horizontal','whisker',1,'PlotStyle','traditional')
    title([probenms{pbnmdex} ' Inst. Displacements Ranked by Fluorescence Intensity (Ascending)'])
    xlabel('Displacements $(nm)$','interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_instdisps_ranked_boxplot'])
    saveas(gcf,[probenms{pbnmdex} '_instdisps_ranked_boxplot.png'])
    
    figure
    boxplot2(sortedinstvelsforboxplot,'labels',sortednamesforboxplot,'orientation','horizontal','whisker',1,'PlotStyle','traditional')
    title([probenms{pbnmdex} ' Inst. Velocities Ranked by Fluorescence Intensity (Ascending)'])
    xlabel('Velocities $(nm s^{-1})$','interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_instvel_ranked_boxplot'])
    saveas(gcf,[probenms{pbnmdex} '_instvel_ranked_boxplot.png'])
    
    figure
    c = categorical(sortednamesforboxplot,sortednamesforboxplot);
    b = 100*sortedloverdispsforboxplot
    barh(c,b)
    xlabel(['Percentage of Tracks over ' num2str(Lover) '$nm$'],'interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_percLover' num2str(Lover) '_ranked_barplot'])
    saveas(gcf,[probenms{pbnmdex} '_percLover' num2str(Lover) '_ranked_barplot.png'])
    
    figure
    c = categorical(sortednamesforboxplot,sortednamesforboxplot);
    b = 100*sortedlovervelsforboxplot
    barh(c,b)
    xlabel(['Percentage of Tracks over ' num2str(Lover) '$nm$'],'interpreter','latex')
    ylabel('Experiment Number and Cell ID')
    set(gca,'FontSize',12)
    savefig([probenms{pbnmdex} '_percVelover' num2str(Lover) '_ranked_barplot'])
    saveas(gcf,[probenms{pbnmdex} '_percVelover' num2str(Lover) '_ranked_barplot.png'])
    
    %fluorescence intensity calc
    celld_f(1,:)= celld(1,:)'./cellp(:,2);
    celld_f(2,:)= celld(2,:)'./cellp(:,2);
    celld_f(3,:)= celld(3,:)'./cellp(:,2);
    celld_f(4,:)= celld(4,:)'./cellp(:,2);
    dispscell_f = dispscell./cellp(:,2);
    velcells_f = velcells./cellp(:,2);
    %other plots
    figure
%     subplot(4,3,1)
    hold on;
    plot(fluornorm(1:10),celld_f(1,1:10),'rx')
    plot(fluornorm(11:20),celld_f(1,11:20),'gx')
    plot(fluornorm(21:30),celld_f(1,21:30),'bx')
    title('fluor vs Med Inst Vel')
    saveas(gcf,[probenms{pbnmdex} 'fluor vs Med Inst Vel'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs Med Inst Vel'])
    figure
%     subplot(4,3,2)
    hold on;
    plot(fluornorm(1:10),celld_f(2,1:10),'rx')
    plot(fluornorm(11:20),celld_f(2,11:20),'gx')
    plot(fluornorm(21:30),celld_f(2,21:30),'bx')
    title('fluor vs med vel')
    saveas(gcf,[probenms{pbnmdex} 'fluor vs med vel'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs med vel'])
    figure
%     subplot(4,3,3)
    hold on;
    plot(fluornorm(1:10),celld_f(3,1:10),'rx')
    plot(fluornorm(11:20),celld_f(3,11:20),'gx')
    plot(fluornorm(21:30),celld_f(3,21:30),'bx')
    title('fluor vs med Inst dis')
    saveas(gcf,[probenms{pbnmdex} 'fluor vs med Inst dis'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs med Inst dis'])
    figure
%     subplot(4,3,4)
    hold on;
    plot(fluornorm(1:10),celld_f(4,1:10),'rx')
    plot(fluornorm(11:20),celld_f(4,11:20),'gx')
    plot(fluornorm(21:30),celld_f(4,21:30),'bx')
    title('fluor vs med dis')
    saveas(gcf,[probenms{pbnmdex} 'fluor vs med dis'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs med dis'])
    figure
%     subplot(4,3,5)
    hold on;
    plot(fluornorm(1:10),dispscell(1:10),'rx')
    plot(fluornorm(11:20),dispscell(11:20),'gx')
    plot(fluornorm(21:30),dispscell(21:30),'bx')
    hold off;
    title(['fluor vs disps over' num2str(Lover)])
    saveas(gcf,[probenms{pbnmdex} 'fluor vs disps over'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs disps over'])
    figure
%     subplot(4,3,6)
    hold on;
    plot(fluornorm(1:10),velcells(1:10),'rx')
    plot(fluornorm(11:20),velcells(11:20),'gx')
    plot(fluornorm(21:30),velcells(21:30),'bx')
    hold off;
    title(['fluor vs vel disp over' num2str(Lover)])
    saveas(gcf,[probenms{pbnmdex} 'fluor vs vel over'],'epsc')
    savefig([probenms{pbnmdex} 'fluor vs vel over'])
    figure
%     subplot(4,3,7)
    hold on;
    plot(cellp(1:10,2),celld(1,1:10),'rx')
    plot(cellp(11:20,2),celld(1,11:20),'gx')
    plot(cellp(21:30,2),celld(1,21:30),'bx')
    title('area vs med Inst vel')
    saveas(gcf,[probenms{pbnmdex} 'area vs med Inst vel'],'epsc')
    savefig([probenms{pbnmdex} 'area vs med Inst vel'])
    figure
%     subplot(4,3,8)
    hold on;
    plot(cellp(1:10,2),celld(2,1:10),'rx')
    plot(cellp(11:20,2),celld(2,11:20),'gx')
    plot(cellp(21:30,2),celld(2,21:30),'bx')
    title('area vs med vel')
    saveas(gcf,[probenms{pbnmdex} 'area vs med vel'],'epsc')
    savefig([probenms{pbnmdex} 'area vs med vel'])
    figure
%     subplot(4,3,9)
    hold on;
    plot(cellp(1:10,2),celld(3,1:10),'rx')
    plot(cellp(11:20,2),celld(3,11:20),'gx')
    plot(cellp(21:30,2),celld(3,21:30),'bx')
    title('area vs med Inst dis')
    saveas(gcf,[probenms{pbnmdex} 'area vs med Inst dis'],'epsc')
    savefig([probenms{pbnmdex} 'area vs med Inst dis'])
    figure
%     subplot(4,3,10)
    hold on;
    plot(cellp(1:10,2),celld(4,1:10),'rx')
    plot(cellp(11:20,2),celld(4,11:20),'gx')
    plot(cellp(21:30,2),celld(4,21:30),'bx')
    title('area vs med dis')
    saveas(gcf,[probenms{pbnmdex} 'area vs med dis'],'epsc')
    savefig([probenms{pbnmdex} 'area vs med dis'])
    figure
%     subplot(4,3,11)
    hold on;
    plot(cellp(1:10,2),dispscell(1:10),'rx')
    plot(cellp(11:20,2),dispscell(11:20),'gx')
    plot(cellp(21:30,2),dispscell(21:30),'bx')
    hold off;
    title(['area vs disps over' num2str(Lover)])
    saveas(gcf,[probenms{pbnmdex} 'area vs disps over'],'epsc')
    savefig([probenms{pbnmdex} 'area vs disps over'])
    figure
%     subplot(4,3,12)
    hold on;
    plot(cellp(1:10,2),velcells(1:10),'rx')
    plot(cellp(11:20,2),velcells(11:20),'gx')
    plot(cellp(21:30,2),velcells(21:30),'bx')
    hold off;
    title(['area vs vel disp over' num2str(Lover)])
%     set(gcf, 'Position',  [100, 100, 800, 800])
%     sgtitle(probenms{pbnmdex})
    saveas(gcf,[probenms{pbnmdex} 'area vs vel over'],'epsc')
    savefig([probenms{pbnmdex} 'area vs vel over'])
%     pause
    cells_meandisps{pbnmdex} = {cells_meandisps;temp_meandisps};
    cells_meanvels{pbnmdex} = {cells_meanvels;temp_meanvels};
    cells_mediandisps{pbnmdex} = {cells_mediandisps;temp_mediandisps};
    cells_medianvels{pbnmdex} = {cells_medianvels;temp_medianvels};
    cells_meaninstdisps{pbnmdex} = {cells_meaninstdisps;temp_meaninstdisps};
    cells_meaninstvels{pbnmdex} = {cells_meaninstvels;temp_meaninstvels};
    cells_medianinstdisps{pbnmdex} = {cells_medianinstdisps;temp_medianinstdisps};
    cells_medianinstvels{pbnmdex} = {cells_medianinstvels;temp_medianinstvels};
    cells_loverdisps{pbnmdex} = {cells_loverdisps;temp_loverdisps};
    cells_lovervels{pbnmdex} = {cells_lovervels;temp_lovervels};
    cells_lovermeandisps{pbnmdex} = {cells_lovermeandisps;temp_lovermeandisps};
    cells_lovermeanvels{pbnmdex} = {cells_lovermeanvels;temp_lovermeanvels};
    cells_lovermediandisps{pbnmdex} = {cells_lovermediandisps;temp_lovermediandisps};
    cells_lovermedianvels{pbnmdex} = {cells_lovermedianvels;temp_lovermedianvels};
end

%% stats
close all;
stats_names = {};
namescount = 1;
stats_meandisps = [];
stats_meanvels = [];
stats_mediandisps = [];
stats_medianvels = [];
stats_meaninstdisps = [];
stats_meaninstvels = [];
stats_medianinstdisps = [];
stats_medianinstvels = [];
stats_loverdisps = [];
stats_lovervels = [];
stats_lovermeandisps = [];
stats_lovermeanvels = [];
stats_lovermediandisps = [];
stats_lovermedianvels = [];
for pdnmdex = 1:length(probenms)
%     disp(probenms{pdnmdex})
    for ndex = 1:length(cells_meandisps{1,pdnmdex}{2,1})
        stats_names{namescount} = probenms{pdnmdex};
        namescount = namescount+1;
    end
    stats_meandisps = [stats_meandisps;cells_meandisps{1,pdnmdex}{2,1}];
    stats_meanvels = [stats_meanvels;cells_meanvels{1,pdnmdex}{2,1}];
    stats_mediandisps = [stats_mediandisps;cells_mediandisps{1,pdnmdex}{2,1}];
    stats_medianvels = [stats_medianvels;cells_medianvels{1,pdnmdex}{2,1}];
    stats_meaninstdisps = [stats_meaninstdisps;cells_meaninstdisps{1,pdnmdex}{2,1}];
    stats_meaninstvels = [stats_meaninstvels;cells_meaninstvels{1,pdnmdex}{2,1}];
    stats_medianinstdisps = [stats_medianinstdisps;cells_medianinstdisps{1,pdnmdex}{2,1}];
    stats_medianinstvels = [stats_medianinstvels;cells_medianinstvels{1,pdnmdex}{2,1}];
    stats_loverdisps = [stats_loverdisps;cells_loverdisps{1,pdnmdex}{2,1}];
    stats_lovervels = [stats_lovervels;cells_lovervels{1,pdnmdex}{2,1}];
    stats_lovermeandisps = [stats_lovermeandisps;cells_lovermeandisps{1,pdnmdex}{2,1}];
    stats_lovermeanvels = [stats_lovermeanvels;cells_lovermeanvels{1,pdnmdex}{2,1}];
    stats_lovermediandisps = [stats_lovermediandisps;cells_lovermediandisps{1,pdnmdex}{2,1}];
    stats_lovermedianvels = [stats_lovermedianvels;cells_lovermedianvels{1,pdnmdex}{2,1}];
end

STATS TESTS

[pmeandisps,~,stats] = anova1(stats_meandisps,stats_names,'off');
[cmeandisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmeandisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmeandisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmeandisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmeandisps(i,1)}...
    ,probenms{cmeandisps(i,2)},cmeandisps(i,3),cmeandisps(i,4),cmeandisps(i,5),cmeandisps(i,6));
end
fclose(fID);


[pmeanvels,~,stats] = anova1(stats_meanvels,stats_names,'off');
[cmeanvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmeanvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmeanvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmeanvels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmeanvels(i,1)}...
    ,probenms{cmeanvels(i,2)},cmeanvels(i,3),cmeanvels(i,4),cmeanvels(i,5),cmeanvels(i,6));
end
fclose(fID);


[pmediandisps,~,stats] = anova1(stats_mediandisps,stats_names,'off');
[cmediandisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmediandisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmediandisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmediandisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmediandisps(i,1)}...
    ,probenms{cmediandisps(i,2)},cmediandisps(i,3),cmediandisps(i,4),cmediandisps(i,5),cmediandisps(i,6));
end
fclose(fID);

[pmedianvels,~,stats] = anova1(stats_medianvels,stats_names,'off');
[cmedianvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmedianvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmedianvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmedianvels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmedianvels(i,1)}...
    ,probenms{cmedianvels(i,2)},cmedianvels(i,3),cmedianvels(i,4),cmedianvels(i,5),cmedianvels(i,6));
end
fclose(fID);


[pmeaninstdisps,~,stats] = anova1(stats_meaninstdisps,stats_names,'off');
[cmeaninstdisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmeaninstdisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmeaninstdisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmeaninstdisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmeaninstdisps(i,1)}...
    ,probenms{cmeaninstdisps(i,2)},cmeaninstdisps(i,3),cmeaninstdisps(i,4),cmeaninstdisps(i,5),cmeaninstdisps(i,6));
end
fclose(fID);

[pmeaninstvels,~,stats] = anova1(stats_meaninstvels,stats_names,'off');
[cmeaninstvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmeaninstvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmeaninstvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmeaninstvels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmeaninstvels(i,1)}...
    ,probenms{cmeaninstvels(i,2)},cmeaninstvels(i,3),cmeaninstvels(i,4),cmeaninstvels(i,5),cmeaninstvels(i,6));
end
fclose(fID);

[pmedianinstdisps,~,stats] = anova1(stats_medianinstdisps,stats_names,'off');
[cmedianinstdisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmedianinstdisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmedianinstdisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmedianinstdisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmedianinstdisps(i,1)}...
    ,probenms{cmedianinstdisps(i,2)},cmedianinstdisps(i,3),cmedianinstdisps(i,4),cmedianinstdisps(i,5),cmedianinstdisps(i,6));
end
fclose(fID);

[pmedianinstvels,~,stats] = anova1(stats_medianinstvels,stats_names,'off');
[cmedianinstvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('pmedianinstvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',pmedianinstvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cmedianinstdisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cmedianinstvels(i,1)}...
    ,probenms{cmedianinstvels(i,2)},cmedianinstvels(i,3),cmedianinstvels(i,4),cmedianinstvels(i,5),cmedianinstvels(i,6));
end
fclose(fID);

[ploverdisps,~,stats] = anova1(stats_loverdisps,stats_names,'off');
[cloverdisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('ploverdisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',ploverdisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(cloverdisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{cloverdisps(i,1)}...
    ,probenms{cloverdisps(i,2)},cloverdisps(i,3),cloverdisps(i,4),cloverdisps(i,5),cloverdisps(i,6));
end
fclose(fID);

[plovervels,~,stats] = anova1(stats_lovervels,stats_names,'off');
[clovervels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('plovervels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',plovervels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(clovervels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{clovervels(i,1)}...
    ,probenms{clovervels(i,2)},clovervels(i,3),clovervels(i,4),clovervels(i,5),clovervels(i,6));
end
fclose(fID);

[plovermeandisps,~,stats] = anova1(stats_lovermeandisps,stats_names,'off');
[clovermeandisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('plovermeandisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',plovermeandisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(clovermeandisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{clovermeandisps(i,1)}...
    ,probenms{clovermeandisps(i,2)},clovermeandisps(i,3),clovermeandisps(i,4),clovermeandisps(i,5),clovermeandisps(i,6));
end
fclose(fID);

[plovermeanvels,~,stats] = anova1(stats_lovermeanvels,stats_names,'off');
[clovermeanvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('plovermeanvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',plovermeanvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(clovermeanvels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{clovermeanvels(i,1)}...
    ,probenms{clovermeanvels(i,2)},clovermeanvels(i,3),clovermeanvels(i,4),clovermeanvels(i,5),clovermeanvels(i,6));
end
fclose(fID);

[plovermediandisps,~,stats] = anova1(stats_lovermediandisps,stats_names,'off');
[clovermediandisps,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('plovermediandisps.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',plovermediandisps)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(clovermediandisps)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{clovermediandisps(i,1)}...
    ,probenms{clovermediandisps(i,2)},clovermediandisps(i,3),clovermediandisps(i,4),clovermediandisps(i,5),clovermediandisps(i,6));
end
fclose(fID);

[plovermedianvels,~,stats] = anova1(stats_lovermedianvels,stats_names,'off');
[clovermedianvels,~,~,gnames] = multcompare(stats,'Alpha',0.05,'CType','tukey-kramer');
fID = fopen('plovermedianvels.csv','w');
fprintf(fID,'pvalue = , %12.5e\n',plovermedianvels)
fprintf(fID,'samp1,samp2,lower95conf,diffmeans,upper95conf,p\n')
for i = 1:length(clovermedianvels)
    fprintf(fID,'%s,%s,%12.5e,%12.5e,%12.5e,%12.5e\n',probenms{clovermedianvels(i,1)}...
    ,probenms{clovermedianvels(i,2)},clovermedianvels(i,3),clovermedianvels(i,4),clovermedianvels(i,5),clovermedianvels(i,6));
end
fclose(fID);

%% save as csv
filename = [foldername 'loverdisps.csv'];
disp(filename)
csvwrite(filename,[cells_loverdisps{1,1}{2,1},cells_loverdisps{1,2}{2,1},cells_loverdisps{1,3}{2,1},cells_loverdisps{1,4}{2,1},cells_loverdisps{1,5}{2,1},cells_loverdisps{1,6}{2,1},cells_loverdisps{1,7}{2,1}]);

filename = [foldername 'lovermeandisps.csv'];
disp(filename)
csvwrite(filename,[cells_lovermeandisps{1,1}{2,1},cells_lovermeandisps{1,2}{2,1},cells_lovermeandisps{1,3}{2,1},cells_lovermeandisps{1,4}{2,1},cells_lovermeandisps{1,5}{2,1},cells_lovermeandisps{1,6}{2,1},cells_lovermeandisps{1,7}{2,1}]);

filename = [foldername 'lovermediandisps.csv'];
disp(filename)
csvwrite(filename,[cells_lovermediandisps{1,1}{2,1},cells_lovermediandisps{1,2}{2,1},cells_lovermediandisps{1,3}{2,1},cells_lovermediandisps{1,4}{2,1},cells_lovermediandisps{1,5}{2,1},cells_lovermediandisps{1,6}{2,1},cells_lovermediandisps{1,7}{2,1}]);

filename = [foldername 'lovermeanvels.csv'];
disp(filename)
csvwrite(filename,[cells_lovermeanvels{1,1}{2,1},cells_lovermeanvels{1,2}{2,1},cells_lovermeanvels{1,3}{2,1},cells_lovermeanvels{1,4}{2,1},cells_lovermeanvels{1,5}{2,1},cells_lovermeanvels{1,6}{2,1},cells_lovermeanvels{1,7}{2,1}]);

filename = [foldername 'lovermedianvels.csv'];
disp(filename)
csvwrite(filename,[cells_lovermedianvels{1,1}{2,1},cells_lovermedianvels{1,2}{2,1},cells_lovermedianvels{1,3}{2,1},cells_lovermedianvels{1,4}{2,1},cells_lovermedianvels{1,5}{2,1},cells_lovermedianvels{1,6}{2,1},cells_lovermedianvels{1,7}{2,1}]);

filename = [foldername 'meandisps.csv'];
disp(filename)
csvwrite(filename,[cells_meandisps{1,1}{2,1},cells_meandisps{1,2}{2,1},cells_meandisps{1,3}{2,1},cells_meandisps{1,4}{2,1},cells_meandisps{1,5}{2,1},cells_meandisps{1,6}{2,1},cells_meandisps{1,7}{2,1}]);

filename = [foldername 'meaninstdisps.csv'];
disp(filename)
csvwrite(filename,[cells_meaninstdisps{1,1}{2,1},cells_meaninstdisps{1,2}{2,1},cells_meaninstdisps{1,3}{2,1},cells_meaninstdisps{1,4}{2,1},cells_meaninstdisps{1,5}{2,1},cells_meaninstdisps{1,6}{2,1},cells_meaninstdisps{1,7}{2,1}]);

filename = [foldername 'medianinstdisps.csv'];
disp(filename)
csvwrite(filename,[cells_medianinstdisps{1,1}{2,1},cells_medianinstdisps{1,2}{2,1},cells_medianinstdisps{1,3}{2,1},cells_medianinstdisps{1,4}{2,1},cells_medianinstdisps{1,5}{2,1},cells_medianinstdisps{1,6}{2,1},cells_medianinstdisps{1,7}{2,1}]);

filename = [foldername 'meaninstvels.csv'];
disp(filename)
csvwrite(filename,[cells_meaninstvels{1,1}{2,1},cells_meaninstvels{1,2}{2,1},cells_meaninstvels{1,3}{2,1},cells_meaninstvels{1,4}{2,1},cells_meaninstvels{1,5}{2,1},cells_meaninstvels{1,6}{2,1},cells_meaninstvels{1,7}{2,1}]);

filename = [foldername 'medianinstvels.csv'];
disp(filename)
csvwrite(filename,[cells_medianinstvels{1,1}{2,1},cells_medianinstvels{1,2}{2,1},cells_medianinstvels{1,3}{2,1},cells_medianinstvels{1,4}{2,1},cells_medianinstvels{1,5}{2,1},cells_medianinstvels{1,6}{2,1},cells_medianinstvels{1,7}{2,1}]);

filename = [foldername 'medianvels.csv'];
disp(filename)
csvwrite(filename,[cells_medianvels{1,1}{2,1},cells_medianvels{1,2}{2,1},cells_medianvels{1,3}{2,1},cells_medianvels{1,4}{2,1},cells_medianvels{1,5}{2,1},cells_medianvels{1,6}{2,1},cells_medianvels{1,7}{2,1}]);

filename = [foldername 'meanvels.csv'];
disp(filename)
csvwrite(filename,[cells_meanvels{1,1}{2,1},cells_meanvels{1,2}{2,1},cells_meanvels{1,3}{2,1},cells_meanvels{1,4}{2,1},cells_meanvels{1,5}{2,1},cells_meanvels{1,6}{2,1},cells_meanvels{1,7}{2,1}]);

filename = [foldername 'mediandisps.csv'];
disp(filename)
csvwrite(filename,[cells_mediandisps{1,1}{2,1},cells_mediandisps{1,2}{2,1},cells_mediandisps{1,3}{2,1},cells_mediandisps{1,4}{2,1},cells_mediandisps{1,5}{2,1},cells_mediandisps{1,6}{2,1},cells_mediandisps{1,7}{2,1}]);






