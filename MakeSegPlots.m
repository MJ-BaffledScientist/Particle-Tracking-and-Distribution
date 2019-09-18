function [BWcombi] = MakeSegPlots(baseNames,anafilenames,expnum,X,xb,yb,polyin,xcenter,ycenter,xs,ys,ratios,matpic)
%function that draws plots and saves them to validate segmentation using
%Otsu and adaptive thresholding. Also, shows boundary plots given data. See
%BoundaryPlots.m function for variable input info.
%Created: 15 Aug 2018 by Daniel S. Han

figure('visible','off')
imshow(X);
hold on;
plot(xb,yb,'r-','LineWidth',3,'DisplayName','Cell Boundary')
plot(polyin,'DisplayName','Cell')
plot(xcenter,ycenter,'g*','DisplayName','Center')
for index = 1:length(ratios)
    plot(xs(:,index),ys(:,index),'-','LineWidth',1.3,'DisplayName',[num2str(ratios(index)*100) '%'])
end
hold off
legend('Location','eastoutside')
saveas(gcf,['cellRatios_' expnum '_' baseNames '_' anafilenames])

% Otsu threshold
% level = graythresh(matpic);
level = multithresh(matpic,2);
BWotsu = imbinarize(matpic,level(2));

% Adaptive threshold
T = adaptthresh(matpic,0.1);
BWadapt = imbinarize(matpic,T);
BWcombi = imbinarize(matpic,T);
%Combined threshold
BWcombi(BWotsu==0)= 0;

figure('visible','off')
subplot(1,3,1)
imshow(BWotsu)
title('Otsu')
subplot(1,3,2)
imshow(BWadapt)
title('Adaptive')
subplot(1,3,3)
imshow(BWcombi)
title('Combination')
saveas(gcf,['blobfinding_' expnum '_' baseNames '_' anafilenames])

figure('visible','off')
subplot(1,2,1)
imshow(BWcombi)
subplot(1,2,2)
hold on;
imshow(X(min(yb):max(yb),min(xb):max(xb),1))
plot(xb-min(xb),yb-min(yb),'r-','LineWidth',1,'DisplayName','Cell Boundary')
xlim([0 max(xb)-min(xb)]);
ylim([0 max(yb)-min(yb)]);
saveas(gcf,['blobpicture_' expnum '_' baseNames '_' anafilenames])