function target = gettarget2DpngANDrandNAU(filename,target,im,fig)
%% Inputs
%   filename
%   target
%       numTargetMax
%       xStep_m
%       yStep_m
%       xOffset_m
%       yOffset_m
%       zOffset_m
%       o_x
%       o_y
%       ampMin
%       ampMax

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
png.xyz_m = [xT,yT,zT]; % xPoint x 3 (x-y-z) x yPoint;
png.xyz_m = reshape(permute(png.xyz_m,[1 3 2]),[],3);

indT = rot90(tMatrix,-1)==true;
png.xyz_m = single(png.xyz_m(indT,:));

png.numTarget = size(png.xyz_m,1);
png.amp = ones(png.numTarget,1);

%% Fit reflectivity image to im domain

% x_m = im.x_m(1):target.xStep_m:im.x_m(end);
% y_m = im.y_m(1):target.yStep_m:im.y_m(end);
% 
% tMatrixPadded = single(flip(tMatrix,1));
% if (length(x_m) > size(tMatrix,2))
%     tMatrixPadded = padarray(tMatrixPadded,[0 floor((length(x_m)-size(tMatrix,2))/2)],0,'pre');
%     tMatrixPadded = padarray(tMatrixPadded,[0 ceil((length(x_m)-size(tMatrix,2))/2)],0,'post');
% else
%     warning("ERROR! Imaging domain size does not include the entire reflectivity function! Increase size!")
% end
% if (length(y_m) > size(tMatrix,1))
%     tMatrixPadded = padarray(tMatrixPadded,[floor((length(y_m)-size(tMatrix,1))/2) 0],0,'pre');
%     tMatrixPadded = padarray(tMatrixPadded,[ceil((length(y_m)-size(tMatrix,1))/2) 0],0,'post');
% else
%     warning("ERROR! Imaging domain size does not include the entire reflectivity function! Increase size!")    
% end
% 
% tMatrixPadded = imresize(tMatrixPadded,[im.numY,im.numX]);
% tMatrixPadded(tMatrixPadded>1) = 1;
% tMatrixPadded(tMatrixPadded<0) = 0;
% 
% target.ideal2D = imgaussfilt(tMatrixPadded,1);

tMatrix_im = interp2(xAxisT,yAxisT,imgaussfilt(single(flip(tMatrix,1)),1.5),im.x_m,im.y_m,'spline',0);
tMatrix_im = imgaussfilt(tMatrix_im,1);
tMatrix_im = tMatrix_im/max(tMatrix_im(:));
tMatrix_im(tMatrix_im>1) = 1;
tMatrix_im(tMatrix_im<0) = 0;

png.ideal2D = tMatrix_im;
png.xyz_m = target.xyz_m;

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
fail = true;

while fail
    target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),target.zOffset_m*ones(target.numTarget,1)]);
    target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(target.numTarget,1);
    R = pdist2(png.xyz_m,target.xyz_m);
    if min(R,[],'all') > 3e-3
        fail = false;
    end
end

%% Create the ideal reflectivity function
target.ideal2D = single(zeros(im.numX,im.numY));
for indTarget = 1:target.numTarget
    temp = single(exp(-(target.o_x)^(-2)*(im.x_m-target.xyz_m(indTarget,1)).^2-(target.o_y)^(-2)*(im.y_m-target.xyz_m(indTarget,2)).^2));
    temp = temp*target.amp(indTarget)/max(temp(:));
    target.ideal2D = target.ideal2D + temp;
end

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D,'FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function");