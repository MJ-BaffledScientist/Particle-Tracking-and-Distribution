function [BWcombi] = mythreshold(I)
%function that thresholds bw image by comparing otsu and adaptive threshold
%and then wiping out pixels in adaptive threshold if they are 0 in the otsu
%also
%Created: 16 Aug 2018 by Daniel S. Han

% Otsu threshold
level = graythresh(I);
BWotsu = imbinarize(I,level);
% Adaptive threshold
T = adaptthresh(I,0.1);
BWadapt = imbinarize(I,T);
BWcombi = imbinarize(I,T);
%Combined threshold
BWcombi(BWotsu==0)= 0;

% figure('visible','on')
% subplot(1,3,1)
% imshow(BWotsu)
% title('Otsu')
% subplot(1,3,2)
% imshow(BWadapt)
% title('Adaptive')
% subplot(1,3,3)
% imshow(BWcombi)
% title('Combination')