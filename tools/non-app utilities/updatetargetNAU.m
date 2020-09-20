function sarData = updatetargetNAU(target,sar,fmcw)
if target.isMIMO
    %% Get distances
    target.R.tx = pdist2(sar.rx.xyz_m,target.xyz_m);
    target.R.rx = pdist2(sar.tx.xyz_m,target.xyz_m);
    target.R.plus = target.R.tx + target.R.rx;
    target.R.times = target.R.tx .* target.R.rx;
    
    %% Get echo signal
    R_T_plus_R_R = target.R.plus;
    isAmplitudeFactor = target.isAmplitudeFactor;
    if isAmplitudeFactor
        amplitudeFactor = target.amp./target.R.times;
    end
    
    k = fmcw.k;
    
    sarData = single(zeros(size(sar.tx.xyz_m,1),fmcw.ADCSamples));
    
    % Create the progress dialog
    d = waitbar(0,"Generating Echo Signal");
    
    for indK = 1:fmcw.ADCSamples
        temp = exp(1j*k(indK)*R_T_plus_R_R);
        if isAmplitudeFactor
            temp = amplitudeFactor .* temp;
        end
        sarData(:,indK) = sum(temp,2);
        
        % Update the progress dialog
        waitbar(indK/length(k),d,"Generating Echo Signal");
    end
    
    delete(d);
    
    % Reshape echo signal
    sarData = reshape(sarData,[sar.size,fmcw.ADCSamples]);
else
    %% Get distances
    target.R = pdist2(sar.vx.xyz_m,target.xyz_m);
    
    %% Get echo signal
    twoR = 2*target.R;
    isAmplitudeFactor = target.isAmplitudeFactor;
    if isAmplitudeFactor
        one_by_R_squared = target.R.^2;
    end
    
    k = fmcw.k;
    
    sarData = single(zeros(size(sar.vx.xyz_m,1),fmcw.ADCSamples));
    
    % Create the progress dialog
    d = waitbar(0,"Generating Echo Signal");
    
    for indK = 1:fmcw.ADCSamples
        temp = exp(1j*k(indK)*twoR);
        if isAmplitudeFactor
            temp = one_by_R_squared .* temp;
        end
        sarData(:,indK) = sum(temp,2);
        
        % Update the progress dialog
        waitbar(indK/length(k),d,"Generating Echo Signal");
    end
    
    delete(d);
    
    % Reshape echo signal
    sarData = reshape(sarData,[sar.size,fmcw.ADCSamples]);
end