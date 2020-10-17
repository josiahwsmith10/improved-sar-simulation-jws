function im = SISO_BPA_2D_array_reconstructImage_2DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
k = single(reshape(fmcw.k,1,1,[]));

%% Get Antenna Positions (primed)
xp_m = reshape(sar.x_m,1,[]);
yp_m = reshape(sar.y_m,[],1);
Z0 = reshape(sar.z_m,1,1,[]);

%% Get Scene Positions (unprimed)
x_m = reshape(im.x_m,1,[]);
y_m = reshape(im.y_m,[],1);

sarImage = zeros(im.numX,im.numY);
%% Use gpuArray if Needed
if target.isGPU
    xp_m = gpuArray(xp_m);
    yp_m = gpuArray(yp_m);
    Z0 = gpuArray(Z0);
    
    x_m = gpuArray(x_m);
    y_m = gpuArray(y_m);
    
    k = gpuArray(k);
    
    sarData = gpuArray(sarData);
    sarImage = gpuArray(sarImage);
end

%% Compute BPA (VERY LONG!)
d = waitbar(0,"Computing BPA, Please Wait!");
count = 0;

numX = im.numX;
numY = im.numY;
isAmplitudeFactor = target.isAmplitudeFactor;
for indX = 1:numX
    for indY = 1:numY
        R = sqrt( (x_m(indX) - xp_m).^2 + (y_m(indY) - yp_m).^2 + Z0.^2);
        if isAmplitudeFactor
            bpaKernel = R.^2 .* exp(-1j*k*2.*R);
        else
            bpaKernel = exp(-1j*k*2.*R);
        end
        sarImage(indX,indY) = sum(sarData.*bpaKernel,'all');
        
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
mesh(h,im.x_m,im.y_m,im.sarImage','FaceColor','interp')
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Reconstructed Image");
view(h,2)
end