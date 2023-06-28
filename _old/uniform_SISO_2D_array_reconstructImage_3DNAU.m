function im = uniform_SISO_2D_array_reconstructImage_3DNAU(sarData,target,fmcw,ant,sar,im,fig)
% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)

%% Compute Wavenumbers
k = single(reshape(fmcw.k,1,1,[]));
L_x = im.nFFTx * sar.xStep_m;
dkX = 2*pi/L_x;
kX = make_kX(dkX,im.nFFTx);

L_y = im.nFFTy * sar.yStep_m;
dkY = 2*pi/L_y;
kY = make_kX(dkY,im.nFFTy)';

kZU = single(reshape(linspace(0,2*max(k),im.nFFTz),1,1,[]));
dkZU = kZU(2) - kZU(1);

if target.isGPU
    reset(gpuDevice)
    k = gpuArray(k);
    kX = gpuArray(kX);
    kY = gpuArray(kY);
    kZU = gpuArray(kZU);
end

kYU = repmat(kY,[1,im.nFFTx,im.nFFTz]);
kXU = repmat(kX,[im.nFFTy,1,im.nFFTz]);
kU = single(1/2 * sqrt(kX.^2 + kY.^2 + kZU.^2));
kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));

%% Compute Focusing Filter
focusingFilter = exp(-1j * kZ * ant.tx.z0_m);
if target.isAmplitudeFactor
    focusingFilter = kZ .* focusingFilter;
end
focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;

%% Zero-Pad Data: s(y,x,k)
if target.isGPU
    sarData = gpuArray(sarData);
end

sarDataPadded = sarData;
sarDataPadded = padarray(sarDataPadded,[floor((im.nFFTy-size(sarData,1))/2) 0],0,'pre');
sarDataPadded = padarray(sarDataPadded,[0 floor((im.nFFTx-size(sarData,2))/2)],0,'pre');

%% Compute FFT across Y & X Dimensions: S(kY,kX,k)
sarDataFFT = fftshift(fftshift(fft(fft(conj(sarDataPadded),im.nFFTy,1),im.nFFTx,2),1),2)/im.nFFTx/im.nFFTy;
clear sarDataPadded sarData

if target.isGPU
    focusingFilter = gpuArray(focusingFilter);
end

%% Stolt Interpolation
% sarImageFFT = zeros(size(kU));
% for ii = 1:size(kU,1)
%     for jj = 1:size(kU,2)
%         tempS = squeeze(sarDataFFT(ii,jj,:) .* focusingFilter(ii,jj,:));
%         tempkU = squeeze(kU(ii,jj,:));
%         sarImageFFT(ii,jj,:) = interpn(k(:),tempS(:),tempkU(:),'linear',0);
%     end
% end
sarImageFFT = interpn(kY(:),kX(:),k(:), sarDataFFT .* focusingFilter ,kYU,kXU,kU,'linear',0);
clear sarDataFFT focusingFilter kY kX k kYU kXU kU

if target.isGPU
    sarImageFFT = gather(sarImageFFT);
    reset(gpuDevice);
    sarIamgeFFT = gpuArray(sarImageFFT);
end

%% Recover Image by IFT: p(y,x,z)
sarImage = single(abs(ifftn(sarImageFFT)));
clear sarDataFFT focusingFilter

%% Reorient Image: p(x,y,z)
sarImage = permute(sarImage,[2,1,3]);

%% Declare Spatial Vectors
x_m = make_x(sar.xStep_m,im.nFFTx);
y_m = make_x(sar.yStep_m,im.nFFTy);
z_m = single(2*pi / (dkZU * im.nFFTz) * (1:im.nFFTz));

%% Resize Image
if max(im.x_m) > max(x_m)
    warning("WARNING: im.nFFTx is too small to see the image!")
end
if max(im.y_m) > max(y_m)
    warning("WARNING: im.nFFTy is too small to see the image!")
end
if max(im.z_m) > max(z_m)
    warning("WARNING: im.nFFTz is too small to see the image!")
end

[X,Y,Z] = ndgrid(im.x_m(:),im.y_m(:),im.z_m(:));
im.sarImage = single(gather(interpn(x_m(:),y_m(:),z_m(:),sarImage,X,Y,Z,'linear',0)));

%% Show Reconstructed Image
if ~fig.isFig
    return;
end
h = fig.SAR2D.h;
f = fig.SAR2D.f;
plotXYZdB(h,f,im.sarImage,im.x_m,im.y_m,im.z_m,[],im.dBMin,"Reconstructed Image",12)

if target.isGPU
    reset(gpuDevice);
end
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