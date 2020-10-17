function im = SISO_BPA_1D_array_reconstructImage_1DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, fmcw.ADCSamples)
k = single(reshape(fmcw.k,1,[]));

%% Get Antenna Positions (primed)
yp_m = reshape(sar.y_m,[],1);
Z0 = reshape(sar.z_m,1,[]);

%% Get Scene Positions (unprimed)
y_m = reshape(im.y_m,[],1);

sarImage = zeros(im.numY,1);
%% Use gpuArray if Needed
if target.isGPU
    yp_m = gpuArray(yp_m);
    Z0 = gpuArray(Z0);
    
    y_m = gpuArray(y_m);
    
    k = gpuArray(k);
    
    sarData = gpuArray(sarData);
    sarImage = gpuArray(sarImage);
end

%% Compute BPA (VERY LONG!)
d = waitbar(0,"Computing BPA, Please Wait!");
count = 0;

numY = im.numY;
isAmplitudeFactor = target.isAmplitudeFactor;
for indY = 1:numY
    R = sqrt((y_m(indY) - yp_m).^2 + Z0.^2);
    if isAmplitudeFactor
        bpaKernel = R.^2 .* exp(-1j*k*2.*R);
    else
        bpaKernel = exp(-1j*k*2.*R);
    end
    sarImage(indY) = sum(sarData.*bpaKernel,'all');
    
    count = count + 1;
    waitbar(count/numel(sarImage),d,"Computing BPA, Please Wait!");
end

delete(d)

im.sarImage = gather(abs(sarImage));

%% Show Reconstructed Image
if ~fig.isFig
    return;
end
h = fig.SAR2D.h;
plot(h,im.y_m,sarImage)
xlabel(h,"y (m)")
ylabel(h,"z (m)")
xlim(h,[im.y_m(1),im.y_m(end)])
title(h,"Reconstructed Image");
end