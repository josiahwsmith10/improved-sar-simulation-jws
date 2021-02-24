function im = SISO_BPA_1D_array_reconstructImage_2DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, fmcw.ADCSamples)
k = single(reshape(fmcw.k,1,[]));

%% Get Antenna Positions (primed)
yp_m = reshape(sar.y_m,[],1);
zp_m = reshape(sar.z_m,1,[]);

%% Get Scene Positions (unprimed)
y_m = reshape(im.y_m,[],1);
z_m = reshape(im.z_m,1,[]);

sarImage = zeros(im.numY,im.numZ);
%% Use gpuArray if Needed
if target.isGPU
    yp_m = gpuArray(yp_m);
    zp_m = gpuArray(zp_m);
    
    y_m = gpuArray(y_m);
    z_m = gpuArray(z_m);
    
    k = gpuArray(k);
    
    sarData = gpuArray(sarData);
    sarImage = gpuArray(sarImage);
end

%% Compute BPA (VERY LONG!)
d = waitbar(0,"Computing BPA, Please Wait!");
count = 0;

numY = im.numY;
numZ = im.numZ;
isAmplitudeFactor = target.isAmplitudeFactor;
for indY = 1:numY
    for indZ = 1:numZ
        R = sqrt((y_m(indY) - yp_m).^2 + (z_m(indZ) - zp_m).^2);
        if isAmplitudeFactor
            bpaKernel = R.^2 .* exp(-1j*k*2.*R);
        else
            bpaKernel = exp(-1j*k*2.*R);
        end
        sarImage(indY,indZ) = sum(sarData.*bpaKernel,'all');
        
        count = count + 1;
        waitbar(count/numel(sarImage),d,"Computing BPA, Please Wait!");
    end
end

delete(d)

im.sarImage = gather(abs(sarImage));

%% Show Reconstructed Image
if ~fig.isFig
    return;
end
h = fig.SAR2D.h;
mesh(h,im.y_m,im.z_m,im.sarImage','FaceColor','interp')
xlabel(h,"y (m)")
ylabel(h,"z (m)")
xlim(h,[im.y_m(1),im.y_m(end)])
ylim(h,[im.z_m(1),im.z_m(end)])
title(h,"Reconstructed Image");
view(h,2)
end