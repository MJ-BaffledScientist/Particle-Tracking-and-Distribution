
%CLEAN SLATE PROTOCOL
clear all;
close all;
warning('off','all');
%loop through all numbers 
expnums = {"Experiment 1","Experiment 2","Experiment 3"};
ncumpixperareaAlltemp = {};
countcells = {0,0,0,0,0,0,0};
disphist = {};
counter = 0;
for edex = 1:length(expnums)
    path = strcat("/Volumes/Mark_Data/Data Files/Organelle Distribution/Vero/Lysosomes/",expnums{edex},"/");
    experiments = {'Untransfected','EGFP-C1','GFP-Tubulin','EMTB-3xGFP','Mn1+2+1+2','Dynein 2219','Kif5c'};
    datasetnames = {};
    for i = 1:length(experiments)
        datasetnames{i,1} = detectfiles(path,[experiments{i} '.mat']);
    end
    ncumpixperareaExp = {};
    ratios = [];
    for i = 1:length(datasetnames)
        disp(strcat(path,datasetnames{i,1}))
        d = load(strcat(path,datasetnames{i,1}));
        %combine all the ncumpixelsPerAreaCell for one experiment
        countcells{i} = countcells{i} + length(d.ncumpixelsPerAreaCell);;
        tempncumpixs = d.ncumpixelsPerAreaCell{1};

        for j = 2:length(d.ncumpixelsPerAreaCell)
            tempncumpixs = tempncumpixs + d.ncumpixelsPerAreaCell{j};
       
        end
        ncumpixperareaExp{i} = tempncumpixs/length(d.ncumpixelsPerAreaCell);
        ncumpixperareaAlltemp{edex,i} = tempncumpixs/length(d.ncumpixelsPerAreaCell);
        ratios = d.ratios;
   
    end
    %plot cdf
    figure
    hold on;
    for i = 1:length(ncumpixperareaExp)
        plot(ratios',ncumpixperareaExp{i},'DisplayName',experiments{i})
    end
    title(expnums{edex})
    xlabel('Percentage Distance to Cell Periphery')
    ylim([0 1])
    ylabel('CDF of Normalised Intensity')
    legend('Location','southeast')
    saveas(gcf,strcat('Vero_CumulativeLysoDist',expnums{edex},'.png'))
    savefig(gcf,strcat('Vero_CumulativeLysoDist',expnums{edex}))
    
end

%% average for all experiments and plot
%calculate averages
ncumpixperareaAll = {};
for j = 1:length(ncumpixperareaAlltemp)
    temp = 0*ncumpixperareaAlltemp{1,1};
    for i = 1:length(expnums)
        temp = temp+ncumpixperareaAlltemp{i,j};
    end
    ncumpixperareaAll{j} = temp/length(expnums);
end

ncumpixperareaAll_std = {};
ncumpixperareaAll_ser = {};
for j = 1:length(ncumpixperareaAlltemp)
    temp = 0*ncumpixperareaAlltemp{1,1};
    for i = 1:length(expnums)
        temp = temp+(ncumpixperareaAlltemp{i,j}-ncumpixperareaAll{j}).^2;
    end
    ncumpixperareaAll_std{j} = sqrt(temp/(length(expnums)-1));
    ncumpixperareaAll_ser{j} = sqrt(temp/(length(expnums)-1))/sqrt(length(expnums));
end


% load p50 control and also plot
pfifty = load('/Volumes/Mark_Data/Data Files/Organelle Distribution/Vero/p50 control/Data_Vero_Lysosomes_Experiment 1_p50.mat');
tempncumpixs_pf = pfifty.ncumpixelsPerAreaCell{1};
for j = 2:length(pfifty.ncumpixelsPerAreaCell)
    tempncumpixs_pf = tempncumpixs_pf + pfifty.ncumpixelsPerAreaCell{j};
end
tempncumpixs_pf = tempncumpixs_pf/length(pfifty.ncumpixelsPerAreaCell);

%% ks test using cdf
Dmatrix = zeros(length(ncumpixperareaAll));
criticalD = zeros(length(ncumpixperareaAll));
significant = zeros(length(ncumpixperareaAll));
for i = 1:length(ncumpixperareaAll)
    for j = 1:length(ncumpixperareaAll)
        Dmatrix(i,j) = max(abs(ncumpixperareaAll{i}-ncumpixperareaAll{j}));
        %for level of significance a = 0.05, c(a) = 1.36. This will change
        %for different a
        criticalD(i,j) = 1.36*sqrt((countcells{i}+countcells{j})/(countcells{i}*countcells{j}));
        if Dmatrix(i,j) > criticalD(i,j)
            significant(i,j) = 1;
        end
    end
end

%plot
figure
hold on;
colors = {[1.00 0.00 0.00],[1.00 0.41 0.16],[1.00 1.00 0.07],[0.00 1.00 0.00],[0.00 0.00 1.00],[0.49 0.18 0.56],[1.00 0.00 1.00]};
for i = 1:length(ncumpixperareaAll)
    e = errorbar(ratios',ncumpixperareaAll{i},ncumpixperareaAll_ser{i},'HandleVisibility','off','Linewidth',0.5,'CapSize',15);
    e.Color = colors{i};
    plot(ratios',ncumpixperareaAll{i},'-','DisplayName',experiments{i},'Linewidth',1.5,'Color',colors{i})
end
% plot(ratios,tempncumpixs_pf,'DisplayName','p50 Control')
xlabel('Percentage Distance to Cell Periphery')
ylim([0 1])
ylabel('CDF of Normalised Intensity')
legend('Location','southeast','FontSize',18)
saveas(gcf,'Vero_CumulativeLysoDist_All.png')
savefig(gcf,'Vero_CumulativeLysoDistAll')

