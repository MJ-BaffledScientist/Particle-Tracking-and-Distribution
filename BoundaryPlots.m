function [xb,yb,polyin,xc,yc,xs,ys,matx,maty,matpic,matx1,maty1,matpic1] = BoundaryPlots(X,X1,ratios)
%This function detects boundary of cell from segmented image and returns
%positions and intensity values.
%Input:
%X is the original image uint8
%X1 is the image with segmented cell only and black outside uint8.
%ratios is an array of how
%Output:
%xb,yb = positions of cell boundary
%polyin = polyshape object of cell
%xc,yc = centroid position from cell boundary
%xs,ys = scaled positions of cell boundary from ratios
%matx,maty = positions of pixels relative to centroid
%matpic = intensity of pixels matching matx,maty
%matx1,maty1 = positions of pixels relative to centroid in list and no
%zeros
%matpic1 = intensity of pixels matching matx,maty in list and no
%zeros
%Created: 15 Aug 2018 by Daniel S. Han
mask = X(:,:,1); % Red
% mask = X(:,:,2); % Green
%find the blank spots
cspace = X1(:,:,1) == 0;
%find the edges
edges = edge(cspace,'sobel');
[y,x] = find(edges);
%find the boundary shape
[k,~] = boundary(x,y);
xb = x(k);
yb = y(k);
%find the center of boundary shape
polyin = polyshape({xb},{yb});
[xc,yc] = centroid(polyin);
%find the various % shrunken shape of cell
xs = ratios.*(xb - xc)+ xc;
ys = ratios.*(yb - yc)+ yc;
% Get the displacements of every pixel and its intensity>>>>>>>>>>>
matpic = im2double(mask(min(yb):max(yb),min(xb):max(xb)));
[matx,maty] = meshgrid(min(xb):max(xb),min(yb):max(yb));
% matdisp = sqrt((maty-yc).^2+(matx-xc).^2);
%NEED TO MAKE ALL THE PIXELS OUTSIDE THE BOUNDARY ZERO SO THAT OTHER
%CELLS DO NOT INTERFERE WITH PLOT
matxlong = matx(:);
matylong = maty(:);
IndexInterior = isinterior(polyin,matxlong,matylong);
Index = reshape(IndexInterior,size(matpic));
matpic = Index.*matpic;
matx = matx - xc;
maty = maty -yc;
% matx = matx.*Index;
% maty = maty.*Index;
%reshape and get rid of 0 intensity values for matpic, matx and
%maty
matpic1 = reshape(matpic,[numel(matpic),1]);
matx1 = reshape(matx,[numel(matx),1]);
maty1 = reshape(maty,[numel(maty),1]);
matx1(matpic==0)=[];
maty1(matpic==0)=[];
matpic1(matpic==0)=[];




