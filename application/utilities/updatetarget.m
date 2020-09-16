function target = updatetarget(app)
%% Get local variables
fmcw = app.fmcw;
sar = app.sar;
ant = app.ant;

%% Test target: PSF
target.numTarget = 1;
target.xyz_m = [0,0,0];
target.xyz_m = single(target.xyz_m);
isAmplitudeFactor = true;

%% Get distances
target.R.tx = pdist2(sar.rx.xyz_m,target.xyz_m);
target.R.rx = pdist2(sar.tx.xyz_m,target.xyz_m);
target.R.plus = target.R.tx + target.R.rx;
target.R.times = target.R.tx .* target.R.rx;

%% Get angles
% sar.tx.AzAngle_rad = zeros(ant.tx.numTx*sar.numY,target.numTarget);
% sar.tx.ElAngle_rad = zeros(ant.tx.numTx*sar.numY,target.numTarget);
% sar.rx.AzAngle_rad = zeros(ant.rx.numRx*sar.numY,target.numTarget);
% sar.rx.ElAngle_rad = zeros(ant.rx.numRx*sar.numY,target.numTarget);
%
% for indTx = 1:ant.vx.numVx*sar.numY
%     for indTarget = 1:target.numTarget
%         x = sar.tx.xyz_m(indTx,1) - target.xyz_m(indTarget,1);
%         y = sar.tx.xyz_m(indTx,2) - target.xyz_m(indTarget,2);
%         z = sar.tx.xyz_m(indTx,3) - target.xyz_m(indTarget,3);
%         sar.tx.AzAngle_rad(indTx,indTarget) = atan2(y,x);
%         sar.tx.ElAngle_rad(indTx,indTarget) = atan2(z,sqrt(x^2+y^2));
%     end
% end
%
% for indRx = 1:ant.vx.numVx*sar.numY
%     for indTarget = 1:target.numTarget
%         x = sar.rx.xyz_m(indRx,1) - target.xyz_m(indTarget,1);
%         y = sar.rx.xyz_m(indRx,2) - target.xyz_m(indTarget,2);
%         z = sar.rx.xyz_m(indRx,3) - target.xyz_m(indTarget,3);
%         sar.rx.AzAngle_rad(indRx,indTarget) = atan2(y,x);
%         sar.rx.ElAngle_rad(indRx,indTarget) = atan2(z,sqrt(x^2+y^2));
%     end
% end

%% Get echo signal
R_T_plus_R_R = target.R.plus;
if isAmplitudeFactor
    one_by_R_T_times_R_R = 1./target.R.times;
end

k = fmcw.k;

clear target
target.sarData = single(zeros(size(sar.tx.xyz_m,1),fmcw.ADCSamples));

% Create the progress dialog
d = uiprogressdlg(app.UIFigure,"Title","Generating Echo Signal",...
    "Message","Please Wait");

for indK = 1:length(k)
    temp = exp(1j*k(indK)*R_T_plus_R_R);
    if isAmplitudeFactor
        temp = one_by_R_T_times_R_R .* temp;
    end
    target.sarData(:,indK) = sum(temp,2);
    
    % Update the progress dialog
    d.Value = indK/length(k);
end

% Reshape echo signal
target.sarData = reshape(target.sarData,[sar.size,length(k)]);

%% Set App variables
app.fmcw = fmcw;
app.ant = ant;
app.sar = sar;

%%
switch app.SARMethodDropDown.Value
    case "Linear"
        
        
    case "Rectilinear"
        
        
    case "Circular"
        
        
    case "Cylindrical"
        
        
end