%code to draw different patterns for virtual controls
%of cell distribution analysis
%Daniel Han 14 Oct 2018
close all;

n = 512;
pic = zeros(n,n);
boundaries = zeros(n,n);
%% square boundaries
s = round(0.1*n);
f = round(0.9*n);
boundaries(s,s:f) = 1;
boundaries(f,s:f) = 1;
boundaries(s:f,s) = 1;
boundaries(s:f,f) = 1;
imwrite(boundaries,'test_square_boundaries.png')
%% circle boundaries
boundaries = zeros(n,n);
x = -n/2:n/2;
y = -n/2:n/2;
[xx yy] = meshgrid(x,y);
u = zeros(size(xx));
u((xx.^2+yy.^2)<200^2)=1;   % radius 200, center at the origin
%find boundaries
B = bwboundaries(u);
figure
plot(B{1}(:,1),B{1}(:,2))
for i = 1:length(B{1})
    boundaries(B{1}(i,1),B{1}(i,2)) = 1;
end
imwrite(boundaries,'test_circle_boundaries.png')
%% pattern draw
%draw a single blob pattern
x = -n/2:n/2;
y = -n/2:n/2;
xoff = 0;
yoff = 0;
r = 10;
xoff2 = 100;
yoff2 = 100;
r2 = 20;
[xx yy] = meshgrid(x,y);
u = zeros(size(xx));
u(((xx+xoff).^2+(yy+yoff).^2)<r^2)=1;   % radius 100, center at the origin
u(((xx+xoff2).^2+(yy+yoff2).^2)<r2^2)=1;
figure
imagesc(u)
imwrite(u,'test_smalltwoblob.png')

%draw a gradient pattern
x = -n/2:n/2;
y = -n/2:n/2;
r = 2;
[xx yy] = meshgrid(x,y);
u = zeros(size(xx));
for i = 0:20:200
    disp(i)
    for theta = 0:360/i:360
        u(((xx+i*cos(pi*theta/180)).^2+(yy+i*sin(pi*theta/180)).^2)<r^2)=1;   % radius 1, center at the origin
    end
end
figure
imagesc(u)
colorbar
imwrite(u,'test_increasingblobstoperiphery.png')

