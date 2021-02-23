function target = gettarget2DpngNAU(filename,target,im,fig)
%% Inputs
%   filename
%   target
%       xStep_m
%       yStep_m
%       xOffset_m
%       yOffset_m
%       zOffset_m

%% Load the iamge in
tMatrix = imread(filename);
tMatrix = tMatrix(:,:,1);
tMatrix(tMatrix>0) = 1;
tMatrix = ~tMatrix;
tMatrix = fliplr(tMatrix);

%% Crete the image domain
[target.sizeY,target.sizeX] = size(tMatrix);
xAxisT = target.xStep_m * (-(target.sizeX-1)/2 : (target.sizeX-1)/2);
yAxisT = target.yStep_m * (-(target.sizeY-1)/2 : (target.sizeY-1)/2);

xAxisT = xAxisT + target.xOffset_m;
yAxisT = yAxisT + target.yOffset_m;
zAxisT = target.zOffset_m;

[zT,xT,yT] = meshgrid(zAxisT,xAxisT,yAxisT);
target.xyz_m = [xT,yT,zT]; % xPoint x 3 (x-y-z) x yPoint;
target.xyz_m = reshape(permute(target.xyz_m,[1 3 2]),[],3);

indT = rot90(tMatrix,-1)==true;
target.xyz_m = single(target.xyz_m(indT,:));

target.numTarget = size(target.xyz_m,1);
target.amp = ones(target.numTarget,1);

%% Fit reflectivity image to im domain
tMatrix_im = interp2(xAxisT,yAxisT,imgaussfilt(single(flip(tMatrix,1)),1.5),im.x_m,im.y_m,'spline',0);
tMatrix_im = imgaussfilt(tMatrix_im,1);
tMatrix_im = tMatrix_im/max(tMatrix_im(:));
tMatrix_im(tMatrix_im>1) = 1;
tMatrix_im(tMatrix_im<0) = 0;

target.ideal2D = tMatrix_im;

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D,'FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function");