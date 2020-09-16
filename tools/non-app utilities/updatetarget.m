function target = updatetarget(target,sar,fmcw)
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

%% Get echo signal
R_T_plus_R_R = target.R.plus;
if isAmplitudeFactor
    one_by_R_T_times_R_R = 1./target.R.times;
end

k = fmcw.k;

clear target
target.sarData = single(zeros(size(sar.tx.xyz_m,1),fmcw.ADCSamples));

% Create the progress dialog
d = waitbar(0,"Generating Echo Signal");

for indK = 1:length(k)
    temp = exp(1j*k(indK)*R_T_plus_R_R);
    if isAmplitudeFactor
        temp = one_by_R_T_times_R_R .* temp;
    end
    target.sarData(:,indK) = sum(temp,2);
    
    % Update the progress dialog
    waitbar(indK/length(k),d,"Generating Echo Signal");
end

delete(d);

% Reshape echo signal
target.sarData = reshape(target.sarData,[sar.size,length(k)]);
return
%%
switch app.SARMethodDropDown.Value
    case "Linear"
        
        
    case "Rectilinear"
        
        
    case "Circular"
        
        
    case "Cylindrical"
        
        
end