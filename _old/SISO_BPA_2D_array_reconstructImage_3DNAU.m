function im = SISO_BPA_2D_array_reconstructImage_3DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
k = single(reshape(fmcw.k,1,1,[]));

%% Get Antenna Positions (primed)
xp_m = reshape(sar.x_m,1,[]);
yp_m = reshape(sar.y_m,[],1);
zp_m = reshape(sar.z_m,1,1,[]);

%% Get Scene Positions (unprimed)
x_m = reshape(im.x_m,1,[]);
y_m = reshape(im.y_m,[],1);
z_m = reshape(im.z_m,1,1,[]);

sarImage = zeros(im.numX,im.numY,im.numZ);
%% Use gpuArray if Needed
if target.isGPU
    xp_m = gpuArray(xp_m);
    yp_m = gpuArray(yp_m);
    zp_m = gpuArray(zp_m);
    
    x_m = gpuArray(x_m);
    y_m = gpuArray(y_m);
    z_m = gpuArray(z_m);
    
    k = gpuArray(k);
    
    sarData = gpuArray(sarData);
    sarImage = gpuArray(sarImage);
end

%% Compute BPA (VERY LONG!)
d = waitbar(0,"Computing BPA, Please Wait!");
count = 0;

numX = im.numX;
numY = im.numY;
numZ = im.numZ;
isAmplitudeFactor = target.isAmplitudeFactor;
for indX = 1:numX
    for indY = 1:numY
        for indZ = 1:numZ
            R = sqrt( (x_m(indX) - xp_m).^2 + (y_m(indY) - yp_m).^2 + (z_m(indZ) - zp_m).^2);
            if isAmplitudeFactor
                bpaKernel = R.^2 .* exp(-1j*k*2.*R);
            else
                bpaKernel = exp(-1j*k*2.*R);
            end
            sarImage(indX,indY,indZ) = sum(sarData.*bpaKernel,'all');
            
            count = count + 1;
            waitbar(count/numel(sarImage),d,"Computing BPA, Please Wait!");
        end
    end
end

delete(d)

im.sarImage = gather(abs(sarImage));

%% Show Reconstructed Image
if ~fig.isFig
    return;
end
h = fig.SAR2D.h;
f = fig.SAR2D.f;
plotXYZdB(h,f,im.sarImage,im.x_m,im.y_m,im.z_m,[],im.dBMin,"BPA Reconstructed Image",12)
end