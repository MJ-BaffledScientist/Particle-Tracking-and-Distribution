function [npix] = pixelspersegment(bwcombi,ratios,xs,ys,xb,yb)
%function that calculates the number of pixels per segment, where segments
%are boundaries of the cell scaled by fractions in ratios.
%bwcombi is the segmented image
%xs,ys are the co-ordinates of the scaled boundaries. Output from
%BoundaryPlots.m
%xb,yb are the co-ordinates of the boundaries. Output from
%BoundaryPlots.m
%Created: 16 Aug 2018 by Daniel S. Han

%storage for the different cell boundary segments
polylevels = cell(length(ratios),1);
npix = nan(length(ratios),1);
%make mesh of picture pixel positions
[meshx,meshy] = meshgrid(min(xb):max(xb),min(yb):max(yb));
%[meshx,meshy] = meshgrid(1:rmask,1:cmask);
meshx = meshx(:);
meshy = meshy(:);
%loop through each ratio value
for pdex = 1:length(ratios)
    polylevels{pdex} = polyshape({xs(:,pdex)},{ys(:,pdex)});
    TFin = isinterior(polylevels{pdex},meshx,meshy);
    [rmask,~] = size(bwcombi);
    ISum = sum(sum(bwcombi((meshx(TFin))*rmask+(meshy(TFin)+1))));
    npix(pdex) = ISum/sum(TFin);
end