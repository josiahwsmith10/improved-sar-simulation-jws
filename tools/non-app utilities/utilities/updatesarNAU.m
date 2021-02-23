function sar = updatesarNAU(sar,ant,fig)
%% Inputs
%   sarIn
%       xStep_m
%       yStep_m
%       thetaMax_deg
%       method
%       numX
%       numY
%       numTheta
%
%   ant
%       tx
%       rx
%       vx

%% Determine SAR Method
switch sar.method
    case "Linear"
        % Verify linearity of MIMO array
        if max(diff([ant.tx.xy(:,1);ant.rx.xy(:,1)])) > 8*eps
            warning("MIMO array must be colinear. Please disable necessary elements.")
            return
        end
        
        sar = getsarAxes(ant,sar);
        
        [sar.X_m,sar.Y_m,sar.Z_m] = ndgrid(sar.x_m,sar.y_m,sar.z_m);
        sar.xyz_m = reshape(cat(3,sar.X_m,sar.Y_m,sar.Z_m),1,[],3);
        sar.xyz_m = repmat(sar.xyz_m,ant.vx.numVx,1,1);
        
        sar.tx.xyz_m = single(sar.xyz_m + ant.tx.xyz_m);
        sar.rx.xyz_m = single(sar.xyz_m + ant.rx.xyz_m);
        sar.vx.xyz_m = single(sar.xyz_m + ant.vx.xyz_m);
        
        sar.tx.xyz_m = reshape(sar.tx.xyz_m,[],3);
        sar.rx.xyz_m = reshape(sar.rx.xyz_m,[],3);
        sar.vx.xyz_m = reshape(sar.vx.xyz_m,[],3);
        
        % Unwrap sar.tx.xyz_m & sar.rx.xyz_m as [numRx,numTx,numY,3]
        sar.size = [ant.rx.numRx,ant.tx.numTx,sar.numY];
        
    case "Rectilinear"
        sar = getsarAxes(ant,sar);
        
        [sar.X_m,sar.Y_m,sar.Z_m] = ndgrid(sar.x_m,sar.y_m,sar.z_m);
        sar.xyz_m = reshape(cat(3,sar.X_m,sar.Y_m,sar.Z_m),1,[],3);
        sar.xyz_m = repmat(sar.xyz_m,ant.vx.numVx,1,1);
        
        sar.tx.xyz_m = single(sar.xyz_m + ant.tx.xyz_m);
        sar.rx.xyz_m = single(sar.xyz_m + ant.rx.xyz_m);
        sar.vx.xyz_m = single(sar.xyz_m + ant.vx.xyz_m);
        
        sar.tx.xyz_m = reshape(sar.tx.xyz_m,[],3);
        sar.rx.xyz_m = reshape(sar.rx.xyz_m,[],3);
        sar.vx.xyz_m = reshape(sar.vx.xyz_m,[],3);
        
        % Unwrap sar.tx.xyz_m & sar.rx.xyz_m as [numRx,numTx,numX,numY,3]
        sar.size = [ant.rx.numRx,ant.tx.numTx,sar.numX,sar.numY];
        
    case "Circular"
        % Verify single element array
        if ant.tx.numTx ~= 1 || ant.rx.numRx ~= 1
            warning("Array must have only 1 Tx and 1 Rx. Please disable necessary elements.")
            return
        end
        
        sar = getsarAxes(ant,sar);
        
        sar.x_m = ant.tx.z0_m*cos(sar.theta_rad);
        sar.y_m = zeros(size(sar.theta_rad));
        sar.z_m = ant.tx.z0_m*sin(sar.theta_rad);
        
        sar.xyz_m = reshape([sar.x_m(:),sar.y_m(:),sar.z_m(:)],1,[],3);
        
        sar.tx.xyz_m = single(sar.xyz_m + ant.tx.xyz_m);
        sar.rx.xyz_m = single(sar.xyz_m + ant.rx.xyz_m);
        sar.vx.xyz_m = single(sar.xyz_m + ant.vx.xyz_m);
        
        sar.tx.xyz_m = reshape(sar.tx.xyz_m,[],3);
        sar.rx.xyz_m = reshape(sar.rx.xyz_m,[],3);
        sar.vx.xyz_m = reshape(sar.vx.xyz_m,[],3);
        
        % Unwrap sar.tx.xyz_m & sar.rx.xyz_m as [numRx,numTx,numTheta,3]
        sar.size = [ant.rx.numRx,ant.tx.numTx,sar.numTheta];
        
    case "Cylindrical"
        % Verify linearity of MIMO array
        if max(diff([ant.tx.xy(:,1);ant.rx.xy(:,1)])) > 8*eps
            warning("MIMO array must be colinear. Please disable necessary elements.")
            return
        end
        
        sar = getsarAxes(ant,sar);
        
        sar.x_m = ant.tx.z0_m*cos(sar.theta_rad);
        sar.z_m = ant.tx.z0_m*sin(sar.theta_rad);
        
        % Use theta as first dimension -> so we obtain s(theta,y')
        sar.X_m = repmat(sar.x_m(:),1,sar.numY);
        sar.Z_m = repmat(sar.z_m(:),1,sar.numY);
        sar.Y_m = repmat(sar.y_m,sar.numTheta,1);
        
        sar.xyz_m = reshape(cat(3,sar.X_m,sar.Y_m,sar.Z_m),1,[],3);
        
        sar.tx.xyz_m = single(sar.xyz_m + ant.tx.xyz_m);
        sar.rx.xyz_m = single(sar.xyz_m + ant.rx.xyz_m);
        sar.vx.xyz_m = single(sar.xyz_m + ant.vx.xyz_m);
        
        sar.tx.xyz_m = reshape(sar.tx.xyz_m,[],3);
        sar.rx.xyz_m = reshape(sar.rx.xyz_m,[],3);
        sar.vx.xyz_m = reshape(sar.vx.xyz_m,[],3);
        
        % Unwrap sar.tx.xyz_m & sar.rx.xyz_m as [numRx,numTx,numTheta,numY,3]
        sar.size = [ant.rx.numRx,ant.tx.numTx,sar.numTheta,sar.numY];
        
end

sar.tx.xyz_m = single(sar.tx.xyz_m);
sar.rx.xyz_m = single(sar.rx.xyz_m);
sar.vx.xyz_m = single(sar.vx.xyz_m);

%% Plot the synthetic aperture
if ~fig.isFig
    return;
end
h = fig.SARAxes.h;
hold(h,'off')
temp = sar.tx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
hold(h,'on')
temp = sar.rx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
xlabel(h,"x (m)")
temp1 = sar.tx.xyz_m(:,1);
temp2 = sar.rx.xyz_m(:,1);
xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
ylabel(h,"z (m)")
temp1 = sar.tx.xyz_m(:,3);
temp2 = sar.rx.xyz_m(:,3);
ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
zlabel(h,"y (m)")
temp1 = sar.tx.xyz_m(:,2);
temp2 = sar.rx.xyz_m(:,2);
zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
title(h,"MIMO Synthetic Aperture")
legend(h,"Tx","Rx")
view(h,3)
daspect(h,[1 1 1])


h = fig.SARVirtualAxes.h;
hold(h,'off')
temp = sar.vx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
xlabel(h,"x (m)")
temp1 = sar.tx.xyz_m(:,1);
temp2 = sar.rx.xyz_m(:,1);
xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
ylabel(h,"z (m)")
temp1 = sar.tx.xyz_m(:,3);
temp2 = sar.rx.xyz_m(:,3);
ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
zlabel(h,"y (m)")
temp1 = sar.tx.xyz_m(:,2);
temp2 = sar.rx.xyz_m(:,2);
zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
title(h,"Virtual Synthetic Aperture")
legend(h,"Vx")
view(h,3)
daspect(h,[1 1 1])
end

function sar = getsarAxes(ant,sar)
%% Update number of steps
sar.xSize_m = sar.numX * sar.xStep_m;
sar.ySize_m = sar.numY * sar.yStep_m;

%% Create synthetic aperture step axes
sar.x_m = (-(sar.numX - 1)/2 : (sar.numX - 1)/2) * sar.xStep_m;
sar.y_m = (-(sar.numY - 1)/2 : (sar.numY - 1)/2) * sar.yStep_m;
sar.theta_rad = linspace(0,sar.thetaMax_deg - sar.thetaMax_deg/sar.numTheta,sar.numTheta)*2*pi/360;
sar.z_m = 0;
end