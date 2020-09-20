function target = gettarget2DrandNAU(target,im,fig)
%% Inputs
%   filename
%   target
%       numTargetMax
%       zOffset_m
%       o_x
%       o_y
%       ampMin
%       ampMax

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),target.zOffset_m*ones(target.numTarget,1)]);
target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(target.numTarget,1);

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
title(h,"Original Reflectivity Function, " + target.numTarget + " targets");
