function sarData = updatetargetNAU(target,sar,fmcw)

isAmplitudeFactor = target.isAmplitudeFactor;
isGPU = target.isGPU;

if target.isMIMO
    %% Get distances
    target.R.tx = pdist2(sar.rx.xyz_m,target.xyz_m);
    target.R.rx = pdist2(sar.tx.xyz_m,target.xyz_m);
    R_T_plus_R_R = target.R.tx + target.R.rx;
    
    %% Amplitude Factor
    if isAmplitudeFactor
        amplitudeFactor = target.amp./(target.R.tx .* target.R.rx);
    else
        amplitudeFactor = 1;
    end
    
else
    %% Get distances
    target.R = pdist2(sar.vx.xyz_m,target.xyz_m);
    
    %% Get echo signal
    R_T_plus_R_R = 2*target.R;
    isAmplitudeFactor = target.isAmplitudeFactor;
    if isAmplitudeFactor
        amplitudeFactor = target.amp./(target.R).^2;
    end
end

if isGPU
    amplitudeFactor = gpuArray(amplitudeFactor);
    R_T_plus_R_R = gpuArray(R_T_plus_R_R);
elseif ~isGPU
    % Create the progress dialog
    d = waitbar(0,"Generating Echo Signal");
end

k = fmcw.k;

sarData = single(zeros(size(sar.tx.xyz_m,1),fmcw.ADCSamples));

numK = fmcw.ADCSamples;
for indK = 1:numK
    temp = exp(1j*k(indK)*R_T_plus_R_R);
    if isAmplitudeFactor
        temp = amplitudeFactor .* temp;
    end
    
    if isGPU
        sarData(:,indK) = single(gather(sum(temp,2)));
    else
        sarData(:,indK) = sum(temp,2);
        % Update the progress dialog
        waitbar(indK/numK,d,"Generating Echo Signal");
    end
end

if ~isGPU
    delete(d);
end

% Reshape echo signal
sarData = reshape(sarData,[sar.size,fmcw.ADCSamples]);
