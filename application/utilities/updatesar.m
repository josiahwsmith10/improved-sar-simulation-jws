function sar = updatesar(app)
%% Get local variables
ant = app.ant;
fmcw = app.fmcw;
fig = app.fig;

%% Update step sizes and theta aperture size
if app.XStepSwitch.Value == "λ"
    sar.xStep_m = app.XStepSizeEditField.Value * fmcw.lambda_m;
elseif app.XStepSwitch.Value == "mm"
    sar.xStep_m = app.XStepSizeEditField.Value * 1e-3;
end

if app.YStepSwitch.Value == "λ"
    sar.yStep_m = app.YStepSizeEditField.Value * fmcw.lambda_m;
elseif app.YStepSwitch.Value == "mm"
    sar.yStep_m = app.YStepSizeEditField.Value * 1e-3;
end

sar.thetaMax_deg = app.ThetaSizedegEditField.Value;

%% Determine SAR Method
switch app.SARMethodDropDown.Value
    case "Linear"
        % Enable necessary fields
        app.XStepSizeEditField.Enable = false;
        app.YStepSizeEditField.Enable = true;
        app.ThetaSizedegEditField.Enable = false;
        app.NumXStepsEditField.Enable = false;
        app.NumYStepsEditField.Enable = true;
        app.NumThetaStepsEditField.Enable = false;
        
        % Verify linearity of MIMO array
        if max(diff([ant.tx.xy(:,1);ant.rx.xy(:,1)])) > 8*eps
            uiconfirm(app.UIFigure,"MIMO array must be colinear. Please disable necessary elements.",'Array Topology Error!',...
                "Options",{'OK'},'Icon','warning');
            app.SARMethodDropDown.Value = "-";
            h = app.SARAxes;
            hold(h,'off')
            scatter3(h,[],[],[])
            h = app.SARVirtualAxes;
            scatter3(h,[],[],[])
            return
        end
        
        % Set no horizontal and no angular scan
        app.NumXStepsEditField.Value = 1;
        app.NumThetaStepsEditField.Value = 1;
        sar = getsarAxes(app,sar);
        
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
        % Enable necessary fields
        app.XStepSizeEditField.Enable = true;
        app.YStepSizeEditField.Enable = true;
        app.ThetaSizedegEditField.Enable = false;
        app.NumXStepsEditField.Enable = true;
        app.NumYStepsEditField.Enable = true;
        app.NumThetaStepsEditField.Enable = false;
        
        % Set no angular scan
        app.NumThetaStepsEditField.Value = 1;
        sar = getsarAxes(app,sar);
        
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
        % Enable necessary fields
        app.XStepSizeEditField.Enable = false;
        app.YStepSizeEditField.Enable = false;
        app.ThetaSizedegEditField.Enable = true;
        app.NumXStepsEditField.Enable = false;
        app.NumYStepsEditField.Enable = false;
        app.NumThetaStepsEditField.Enable = true;
        
        % Verify single element array
        if ant.tx.numTx ~= 1 || ant.rx.numRx ~= 1
            uiconfirm(app.UIFigure,"Array must have only 1 Tx and 1 Rx. Please disable necessary elements.",'Array Topology Error!',...
                "Options",{'OK'},'Icon','warning');
            app.SARMethodDropDown.Value = "-";
            h = app.SARAxes;
            hold(h,'off')
            scatter3(h,[],[],[])
            h = app.SARVirtualAxes;
            scatter3(h,[],[],[])
            return
        end
        
        % Set no horizontal and no vertical scan
        app.NumXStepsEditField.Value = 1;
        app.NumYStepsEditField.Value = 1;
        sar = getsarAxes(app,sar);
        
        sar.x_m = app.TransmitterZmEditField.Value*cos(sar.theta_rad);
        sar.y_m = zeros(size(sar.theta_rad));
        sar.z_m = app.TransmitterZmEditField.Value*sin(sar.theta_rad);
        
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
        % Enable necessary fields
        app.XStepSizeEditField.Enable = false;
        app.YStepSizeEditField.Enable = true;
        app.ThetaSizedegEditField.Enable = true;
        app.NumXStepsEditField.Enable = false;
        app.NumYStepsEditField.Enable = true;
        app.NumThetaStepsEditField.Enable = true;
        
        % Verify linearity of MIMO array
        if max(diff([ant.tx.xy(:,1);ant.rx.xy(:,1)])) > 8*eps
            uiconfirm(app.UIFigure,"MIMO array must be colinear. Please disable necessary elements.",'Array Topology Error!',...
                "Options",{'OK'},'Icon','warning');
            app.SARMethodDropDown.Value = "-";
            h = app.SARAxes;
            hold(h,'off')
            scatter3(h,[],[],[])
            h = app.SARVirtualAxes;
            scatter3(h,[],[],[])
            return
        end
        
        % Set no horizontal scan
        app.NumXStepsEditField.Value = 1;
        sar = getsarAxes(app,sar);
        sar.x_m = app.TransmitterZmEditField.Value*cos(sar.theta_rad);
        sar.z_m = app.TransmitterZmEditField.Value*sin(sar.theta_rad);
        
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

%% Set SAR method
sar.method = app.SARMethodDropDown.Value;

if sar.method == "-"
    h = app.SARAxes;
    hold(h,'off')
    scatter3(h,[],[],[])
    h = app.SARVirtualAxes;
    scatter3(h,[],[],[])
    return
end

%% Plot the synthetic aperture
h = app.SARAxes;
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

h = app.SARVirtualAxes;
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

%% Set App variables
app.fig = fig;
end

function sar = getsarAxes(app,sar)
%% Update number of steps
sar.numX = app.NumXStepsEditField.Value;
sar.numY = app.NumYStepsEditField.Value;
sar.numTheta = app.NumThetaStepsEditField.Value;

app.XSizemEditField.Value = sar.numX * sar.xStep_m;
app.YSizemEditField.Value = sar.numY * sar.yStep_m;

%% Create synthetic aperture step axes
sar.x_m = (-(sar.numX - 1)/2 : (sar.numX - 1)/2) * sar.xStep_m;
sar.y_m = (-(sar.numY - 1)/2 : (sar.numY - 1)/2) * sar.yStep_m;
sar.theta_rad = linspace(0,sar.thetaMax_deg - sar.thetaMax_deg/sar.numTheta,sar.numTheta)*2*pi/360;
sar.z_m = app.TransmitterZmEditField.Value;
end