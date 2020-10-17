function im = uniform_SISO_RMA_1D_array_reconstructImage_2DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, fmcw.ADCSamples)

%% Compute Wavenumbers 
k = single(reshape(fmcw.k,1,[]));

L_y = im.nFFTy * sar.yStep_m;
dkY = 2*pi/L_y;
kY = make_kX(dkY,im.nFFTy)';

kZU = single(reshape(linspace(0,2*max(k),im.nFFTz),1,[]));
dkZU = kZU(2) - kZU(1);

if target.isGPU
    k = gpuArray(k);
    kY = gpuArray(kY);
    kZU = gpuArray(kZU);
end

kYU = repmat(kY,[1,im.nFFTz]);
kU = single(1/2 * sqrt(kY.^2 + kZU.^2));

if target.isGPU
    sarData = gpuArray(sarData);
end

sarDataPadded = sarData;
sarDataPadded = padarray(sarDataPadded,[floor((im.nFFTy-size(sarData,1))/2) 0],0,'pre');

sarDataFFT = fftshift(fft(conj(sarDataPadded),im.nFFTy,1),1)/im.nFFTy;
clear sarDataUpsample sarData

sarImageFFT = interpn(kY,k,sarDataFFT,kYU,kU,'linear',0);
clear sarDataFFT

sarImage = single(abs(ifft2(sarImageFFT)));
clear sarImageFFT

y_m = make_x(sar.yStep_m,im.nFFTy);
z_m = single(2*pi / (dkZU * im.nFFTz) * (1:im.nFFTz));

% Resize Image
if max(im.y_m) > max(y_m)
    warning("WARNING: im.nFFTy is too small to see the image!")
end
if max(im.z_m) > max(z_m)
    warning("WARNING: im.nFFTz is too small to see the image!")
end

im.sarImage = interpn(y_m,z_m,single(gather(sarImage)),im.y_m,im.z_m,'spline',0);

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

function x = make_x(xStep_m,nFFTx)
x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
x = single(x);
end

function kX = make_kX(dkX,nFFTx)
if mod(nFFTx,2)==0
    kX = dkX * ( -nFFTx/2 : nFFTx/2-1 );
else
    kX = dkX * ( -(nFFTx-1)/2 : (nFFTx-1)/2 );
end
kX = single(kX);
end