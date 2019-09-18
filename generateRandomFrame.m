clear all;
close all;
%% create simulation of square image of size 'nsize'
nsize = [1040 1392];
%uniform random numbers
I = rand(nsize);
%poisson random numbers
Ipos = poissrnd(1,nsize);
%cluster of gaussians
Ig = zeros(nsize);
x = [1:nsize(1)]';
y = [1:nsize(2)]';
yc = 223;
xc = 403;
[X,Y] = meshgrid(x,y);
mu = [[xc yc];[xc+400 yc]];
F = mvnpdf([X(:) Y(:)],[xc+200,yc],[300 300]);
F = reshape(F,length(x),length(y));
Ig = Ig+255*F;

%completely uniform pattern
Isl = zeros(nsize);
Isl(1:3:nsize(1),1:3:nsize(2)) = 1;

%% Plot patterns
figure
imshow(Ig)

imwrite(Ig,'testcentercluster.tif')



