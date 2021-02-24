function im = uniform_SISO_FFT_2D_array_reconstructImage_2DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)

%% Compute Wavenumbers 
k = reshape(fmcw.k,1,1,[]);
L_x = im.nFFTx * sar.xStep_m;
dkX = 2*pi/L_x;
kX = make_kX(dkX,im.nFFTx);

L_y = im.nFFTy * sar.yStep_m;
dkY = 2*pi/L_y;
kY = make_kX(dkY,im.nFFTy)';

if target.isGPU
    k = gpuArray(k);
    kX = gpuArray(kX);
    kY = gpuArray(kY);
end

kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));

%% Compute Focusing Filter
focusingFilter = kZ .* exp(-1j * kZ * abs(ant.tx.z0_m - target.zOffset_m));
focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;
clear k kX kY

if target.isGPU
    sarData = gpuArray(sarData);
end

sarDataPadded = sarData;
sarDataPadded = padarray(sarDataPadded,[floor((im.nFFTy-size(sarData,1))/2) 0],0,'pre');
sarDataPadded = padarray(sarDataPadded,[0 floor((im.nFFTx-size(sarData,2))/2)],0,'pre');

sarDataFFT = fftshift(fftshift(fft(fft(sarDataPadded,im.nFFTy,1),im.nFFTx,2),1),2)/im.nFFTx/im.nFFTy;
clear sarDataUpsample sarData

if target.isGPU
    focusingFilter = gpuArray(focusingFilter);
end

sarImage = single(abs(sum(ifft(ifft(sarDataFFT .* focusingFilter,[],1),[],2),3)))';
clear sarDataFFT focusingFilter

% sarImage = flip(sarImage,2);

x_m = make_x(sar.xStep_m,im.nFFTx);
y_m = make_x(sar.yStep_m,im.nFFTy);

% Resize Image
if max(im.x_m) > max(x_m)
    warning("WARNING: im.nFFTx is too small to see the image!")
end
if max(im.y_m) > max(y_m)
    warning("WARNING: im.nFFTy is too small to see the image!")
end
im.sarImage = interpn(x_m,y_m,single(gather(sarImage)),im.x_m,im.y_m,'spline',0);

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