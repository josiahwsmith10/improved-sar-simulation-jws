function target = gettargetXYZpngANDrandNAU(filename,target,im,fig)
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
%       o_z
%       ampMin
%       ampMax
%       ampAdjust
%       downSample

png = getpngXYZ(filename,target);
target = gettargetXYZ(target,png,im);

%% Create the ideal reflectivity function
target.ideal3D = single(zeros(im.numX,im.numY,im.numZ));

x_m = im.x_m;
y_m = im.y_m;
z_m = im.z_m;
o_x = target.o_x;
o_y = target.o_y;
o_z = target.o_z;
xyz_m = cat(1,target.xyz_m,png.xyz_m);
amp = cat(2,target.amp,ones(size(png.amp)));
target.numTarget = target.numTarget + png.numTarget;

if target.isGPU
    x_m = gpuArray(x_m);
    y_m = gpuArray(y_m);
    z_m = gpuArray(z_m);
    o_x = gpuArray(o_x);
    o_y = gpuArray(o_y);
    o_z = gpuArray(o_z);
    xyz_m = gpuArray(xyz_m);
    amp = gpuArray(amp);
end

for indTarget = 1:target.numTarget
    temp = single(exp(-(o_x)^(-2)*(x_m-xyz_m(indTarget,1)).^2 -(o_y)^(-2)*(y_m-xyz_m(indTarget,2)).^2 -(o_z)^(-2)*(z_m-xyz_m(indTarget,3)).^2));
    tempMax = max(temp(:));
    if tempMax == 0
        tempMax = 1;
    end
    
    temp = temp*amp(indTarget)/tempMax;
    target.ideal3D = target.ideal3D + temp;
end
if abs(max(target.ideal3D(:))-max(target.amp)) > 5e-2
    warning("Amplitude error! The ideal image is distorted!")
end
target.ideal3D(target.ideal3D>1) = 1;
target.ideal3D(target.ideal3D<0) = 0;

target.ideal3D = single(gather(target.ideal3D));

target.xyz_m = cat(1,target.xyz_m,png.xyz_m);
target.amp = cat(2,target.amp,png.amp);
%% Show the reflectivity function
if ~fig.isFig
    return;
end
h = fig.Target2D.h;
f = fig.Target2D.f;
plotXYZdB(h,f,target.ideal3D,im.x_m,im.y_m,im.z_m,[],target.dBMin,"Original Reflectivity Function, " + target.numTarget + " targets",12)
end

function png = getpngXYZ(filename,target)
if isempty(filename)
    png.numTarget = 0;
    png.xyz_m = [];
    png.amp = [];
    return;
end

%% Load the image in
tMatrix = imread(filename);
tMatrix = tMatrix(:,:,1);
tMatrix(tMatrix<64) = 0;
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

png.xyz_m = downsample(png.xyz_m,target.downSample);

png.numTarget = size(png.xyz_m,1);
png.amp = ones(1,png.numTarget)*target.ampAdjust;
end

function target = gettargetXYZ(target,png,im)
if target.numTargetMax == 0
    target.numTarget = 0;
    target.xyz_m = [];
    target.amp = [];
    return;
end
%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.y_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),im.z_m(1) + (im.z_m(end)-im.z_m(1))*rand(target.numTarget,1)]);
target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(1,target.numTarget);

fail = true;
tic

while fail
    if isempty(png.xyz_m)
        R = 1e3;
    else
        R = pdist2(png.xyz_m,target.xyz_m);
    end
    R_same = pdist2(target.xyz_m,target.xyz_m) + 1e3*eye(target.numTarget);
    indGood = min(R,[],1)>(target.o_x*4) & min(R_same,[],1)>(target.o_x*4);
    numGood = sum(indGood);
    if numGood == target.numTarget
        fail = false;
    else
        xyz_m_good = target.xyz_m(indGood,:);
        xyz_m_new = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),im.y_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),im.z_m(1) + (im.z_m(end)-im.z_m(1))*rand(target.numTarget-size(xyz_m_good,1),1)]);
        target.xyz_m = cat(1,xyz_m_good,xyz_m_new);
        
        if toc > 10
            indGood = min(R,[],1)>(target.o_x*4);
            target.xyz_m = target.xyz_m(indGood,:);
            target.amp = target.amp(indGood);
            target.numTarget = sum(indGood);
            warning("Could not place targets correctly in 10s, reducing number of targets")
            break;
        end
    end
end
end