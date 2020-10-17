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
tMatrix_im = interp2(xAxisT,yAxisT,imgaussfilt(single(flip(tMatrix,1)),1.5),im.x_m,im.y_m,'spline',0);
tMatrix_im = imgaussfilt(tMatrix_im,1);
tMatrix_im = tMatrix_im/max(tMatrix_im(:));
tMatrix_im(tMatrix_im>1) = 1;
tMatrix_im(tMatrix_im<0) = 0;

png.ideal2D = tMatrix_im;

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
fail = true;

tic
target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),target.zOffset_m*ones(target.numTarget,1)]);
while fail
    target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(target.numTarget,1);
    R = pdist2(png.xyz_m,target.xyz_m);
    if min(R,[],'all') > 3e-3
        fail = false;
    else
        indGood = min(R,[],1)>3e-3;
        numGood = sum(indGood);
        xyz_m_good = target.xyz_m(indGood,:);
        xyz_m_new = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),target.zOffset_m*ones(target.numTarget-size(xyz_m_good,1),1)]);
        target.xyz_m = cat(1,xyz_m_good,xyz_m_new);
%         disp(target.numTarget + " total," + (target.numTarget - numGood) + " to go")
        
        if toc > 10
            indGood = min(R,[],1)>3e-3;
            target.xyz_m = target.xyz_m(indGood,:);
            target.amp = target.amp(indGood);
            target.numTarget = sum(indGood);
            return;
        end
    end
end

%% Create the ideal reflectivity function
target.ideal2D = single(zeros(im.numX,im.numY));
for indTarget = 1:target.numTarget
    temp = single(exp(-(target.o_x)^(-2)*(im.x_m-target.xyz_m(indTarget,1)).^2-(target.o_y)^(-2)*(im.y_m-target.xyz_m(indTarget,2)).^2));
    temp = temp*target.amp(indTarget)/max(temp(:));
    target.ideal2D = target.ideal2D + temp;
end
target.ideal2D(target.ideal2D>1) = 1;
target.ideal2D(target.ideal2D<0) = 0;

target.ideal2D = target.ideal2D + png.ideal2D;
target.xyz_m = cat(1,target.xyz_m,png.xyz_m);
target.numTarget = target.numTarget + png.numTarget;

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D,'FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function, " + target.numTarget + " targets");