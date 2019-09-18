%Code for determining cell boundaries and distributions of transfected
%organelles.
%7 Aug 2018
%Daniel Han
%include all hot dots and positions and boundaries in the data output.
%%
%CLEAN SLATE PROTOCOL
clear all;
warning('off','all');
%loop through all numbers 
expnums = {'Experiment 1'};
path = strcat('/Users/mqbprmj2/Documents/Imaging/Organelle Distribution/Vero Lysosome Distribution/p50/',expnums{1},'/');
experiments = {'p50'};
datasetnames = {};
% for i = 1:length(experiments)
%     datasetnames{i,1} = detectfiles(path,[experiments{i} '.mat']);
% end
% ncumpixperareaExp = {};
% ratios = [];
% for i = 1:length(datasetnames)
%     disp(strcat(path,datasetnames{i,1}))
%     d = load(strcat(path,datasetnames{i,1}));
%     %combine all the ncumpixelsPerAreaCell for one experiment
%     tempncumpixs = d.ncumpixelsPerAreaCell{1};
%     for j = 2:length(d.ncumpixelsPerAreaCell)
%         tempncumpixs = tempncumpixs + d.ncumpixelsPerAreaCell{j};
%     end
%     ncumpixperareaExp{i} = tempncumpixs/length(d.ncumpixelsPerAreaCell);
%     ratios = d.ratios;
% %     ncumpixperareaExp{i} = d.ncumpixelsPerAreaCell;
% end
% 
% figure
% hold on;
% for i = 1:length(ncumpixperareaExp)
%     plot(ratios,ncumpixperareaExp{i},'DisplayName',experiments{i})
% end
% legend()
% expnums = {'Experiment 3'};
%loop through all the experiments

% experiments = {'EMTB-3xGFP','GFP-Tubulin'};
%cell with elements to remove from number searching for each color combine
%file
isname = {'p50'};
% isname = {'EMTB-3xGFP','GFP-Tubulin'};
% for expnames = 1
% for expnum = 1:length(expnums)
%     experimentnumber = expnums{expnum};
%     for expnames = 1:length(experiments)
%         % get images
%         get the full image of Color Combine
%         basefolder = ['/Volumes/Mark Data/Images/Organelle Distribution/HeLaM/Lysosomes/' experimentnumber '/'];
%         basefolder = 'C:\Users\MBCX3DH2\Dropbox (The University of Manchester)\Organelle Distribution\HeLaM Lysosome-Golgi\';
%         type = experiments{expnames};
%         disp(experiments{expnames});
%         folder = [basefolder type '\Lysosomes Composites\'];
%         folder = [basefolder type '/Colour Combine Lysosomes/'];
%         anafolder = [basefolder 'Analysis/' type '/'];
%         baseNames = detectfiles(folder,'.tif');
%         for basedex = 1:length(baseNames)
%             filename = fullfile(folder, baseNames{basedex});
%             Read in the image for Color Combine
%             [X,map] = imread(filename);
%             get the segmented folders
%             disp(anafolder);
%             disp(type(isstrprop(type,'digit')))
%             findname = erase(baseNames{basedex},isname{expnames});
%             disp(findname)
%             disp(findname(isstrprop(findname,'digit')))
%             anadirname = [anafolder 'Image ' findname(isstrprop(findname,'digit')) '/'];
%             disp(anadirname)
%             anafilenames = detectfiles(anadirname,'.jpg');
%                 for anadex = 1
%         end
%     end
% end
% pause
for expnum = 1:length(expnums)
    experimentnumber = expnums{expnum};
    for expnames = 1:length(experiments)
        %make cell arrays for all the variables that need to be stored
        ccounter = 0;
        displacementsCell = {};
        totalon = {};
        npixelsPerAreaCell = {};
        ncumpixelsPerAreaCell = {};
        ripleyKdata ={};
        %% get images
        %get the full image of Color Combine
        basefolder = ['/Users/mqbprmj2/Documents/Imaging/Organelle Distribution/Vero Lysosome Distribution/p50/' experimentnumber '/'];
    %     basefolder = 'C:\Users\MBCX3DH2\Dropbox (The University of Manchester)\Organelle Distribution\HeLaM Lysosome-Golgi\';
        type = experiments{expnames};
        disp(experiments{expnames});
    %     folder = [basefolder type '\Lysosomes Composites\'];
        folder = [basefolder type '/Colour Combine Lysosomes/'];
        anafolder = [basefolder 'Analysis/' type '/'];
        baseNames = detectfiles(folder,'.tif');
        % for basedex = 1
        for basedex = 1:length(baseNames)
            filename = fullfile(folder, baseNames{basedex});
            %Read in the image for Color Combine
            [X,map] = imread(filename);
            %get the segmented folders
            findname = erase(baseNames{basedex},isname{expnames});
            disp(anafolder);
            disp(type(isstrprop(type,'digit')))
            disp(findname(isstrprop(findname,'digit')))
            anadirname = [anafolder 'Image ' findname(isstrprop(findname,'digit')) '/'];
            disp(anadirname)
            anafilenames = detectfiles(anadirname,'.jpg');
            %     for anadex = 1
            for anadex = 1:length(anafilenames)
                disp(anafilenames{anadex})
                %Read in the images for Segmented Cells
                [X1,map1] = imread([anadirname anafilenames{anadex}]);
                %% Begin Analysis
                %ratios = 0.5;
                ratios = 0.1:0.05:1;
                [xb,yb,polyin,xcenter,ycenter,xs,ys,matx,maty,matpic,matx1,maty1,matpic1] = BoundaryPlots(X,X1,ratios);
                %% Do the Segmentation
                close all;
                bwcombi = MakeSegPlots(baseNames{basedex}(1:end-4),anafilenames{anadex},experimentnumber,X,xb,yb,polyin,xcenter,ycenter,xs,ys,ratios,matpic);
                %shift centers so that it is correct for bwcombi
                xcenter = xcenter-min(xb);
                ycenter = ycenter-min(yb);
                xs = xs - min(xb);
                ys = ys - min(yb);
                xb = xb - min(xb);
                yb = yb - min(yb);
                %% Measure the number of 'on' pixels
                npix = pixelspersegment(bwcombi,ratios,xs,ys,xb,yb);
                npixcum = cumpixelspersegment(bwcombi,ratios,xs,ys,xb,yb);
                figure('visible','off')
                plot(100*ratios,npix,'r.-');
                xlabel('Percentage of cell boundary from centroid (\%)','interpreter','latex')
                ylabel('Ratio of pixels per area','interpreter','latex')
                saveas(gcf,['pixelsPerSegPlot_' experimentnumber '_' baseNames{basedex}(1:end-4) '_' anafilenames{anadex}])
                %plot cumulative distribution
                figure('visible','off')
                plot(100*ratios,npixcum,'r.-');
                xlabel('Percentage of cell boundary from centroid (\%)','interpreter','latex')
                ylabel('Cumulative Ratio of pixels per area','interpreter','latex')
                %         ylim([0.58 0.64]) % for uniform
                %         ylim([0.07 0.13]) % for square lattice
                saveas(gcf,['C_cumpixelsPerSegPlot_' experimentnumber '_' baseNames{basedex}(1:end-4) '_' anafilenames{anadex}])
                %% Collect all displacements from 'on' pixels
%                 [displ,xK,K] = pixelscount(bwcombi,xcenter,ycenter);
%                 figure('visible','off')
%                 histogram(displ)
%                 xlabel('Distance from centroid (pix.)','interpreter','latex')
%                 ylabel('Counts','interpreter','latex')
%                 saveas(gcf,['DispHistPlot_' experimentnumber '_' baseNames{basedex}(1:end-4) '_' anafilenames{anadex}])
                %% Store the analysis data
                ccounter = ccounter+1;
%                 displacementsCell{ccounter} = displ;
                npixelsPerAreaCell{ccounter} = npix;
                ncumpixelsPerAreaCell{ccounter} = npixcum;
                [yn,xn] = find(bwcombi);
                totalon{ccounter} = length(xn);
%                 ripleyKdata{ccounter,1} = xK;
%                 ripleyKdata{ccounter,2} = K;
            end
        end
        save(['Data_Vero_Lysosomes_' experimentnumber '_' type]);
    end
end