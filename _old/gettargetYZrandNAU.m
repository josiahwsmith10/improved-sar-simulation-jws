function target = gettargetYZrandNAU(target,im,fig)
%% Inputs
%   filename
%   target
%       numTargetMax
%       zOffset_m
%       o_y
%       o_z
%       ampMin
%       ampMax

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
target.xyz_m = single([zeros(target.numTarget,1),im.y_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),im.z_m(1) + (im.z_m(end)-im.z_m(1))*rand(target.numTarget,1)]);
target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(1,target.numTarget);

fail = true;
tic;

while fail
    R_same = pdist2(target.xyz_m,target.xyz_m) + 1e3*eye(target.numTarget);
    indGood = min(R_same,[],1)>(target.o_y*4);
    numGood = sum(indGood);
    if numGood == target.numTarget
        fail = false;
    else
        xyz_m_good = target.xyz_m(indGood,:);
        xyz_m_new = single([zeros(target.numTarget-size(xyz_m_good,1),1),im.y_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),im.z_m(1) + (im.z_m(end)-im.z_m(1))*rand(target.numTarget-size(xyz_m_good,1),1)]);
        target.xyz_m = cat(1,xyz_m_good,xyz_m_new);
        
        if toc > 10
            indGood = min(R_same,[],1)>(target.o_y*4);
            target.xyz_m = target.xyz_m(indGood,:);
            target.amp = target.amp(indGood);
            target.numTarget = sum(indGood);
            warning("Could not place targets correctly in 10s, reducing number of targets")
            break;
        end
    end
end

%% Create the ideal reflectivity function
target.ideal2D = single(zeros(im.numY,im.numZ));
for indTarget = 1:target.numTarget
    temp = single(exp(-(target.o_y)^(-2)*(im.y_m-target.xyz_m(indTarget,2)).^2-(target.o_z)^(-2)*(im.z_m-target.xyz_m(indTarget,3)).^2));
    tempMax = max(temp(:));
    if tempMax == 0
        tempMax = 1;
    end
    
    temp = temp*target.amp(indTarget)/tempMax;
    target.ideal2D = target.ideal2D + temp;
    
    if sum(isnan(target.ideal2D))
        warning("oops")
    end
end
if abs(max(target.ideal2D(:))-max(target.amp)) > 3e-2
    warning("Amplitude error! The ideal image is distorted!")
end
target.ideal2D(target.ideal2D>1) = 1;
target.ideal2D(target.ideal2D<0) = 0;

%% Show the reflectivity function
if ~fig.isFig
    return;
end
h = fig.Target2D.h;
mesh(h,im.y_m,im.z_m,target.ideal2D','FaceColor','interp')
view(h,2)
xlabel(h,"y (m)")
ylabel(h,"z (m)")
xlim(h,[im.y_m(1),im.y_m(end)])
ylim(h,[im.z_m(1),im.z_m(end)])
title(h,"Original Reflectivity Function, " + target.numTarget + " targets");
