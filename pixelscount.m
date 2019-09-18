function [disp,xK,K] = pixelscount(bwcombi,xcenter,ycenter)
[r,~] = size(bwcombi);
[y,x] = find(bwcombi);
y = r-y;
xc = xcenter;
yc = r-ycenter;

disp = sqrt((x-xc).^2+(y-yc).^2);

%% calculate Ripley K
dataXY = [x,y];
xK = 0:0.1*max(disp):max(disp);
K = ripleyk(dataXY,xK,[min(x), max(x), min(y), max(y)],0);

% figure 
% subplot(1,2,1)
% imshow(bwcombi);
% subplot(1,2,2)
% hold on;
% plot(x,y,'g.');
% plot(xcenter,r-ycenter,'r*')
