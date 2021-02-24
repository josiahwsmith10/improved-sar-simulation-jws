function im = uniform_SISO_2D_array_reconstructImage_2DNAU(sarData,target,fmcw,ant,sar,im)
% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)

k = reshape(fmcw.k,1,1,[]);

sarDataUpsample = sarData;
sarDataUpsample = padarray(sarDataUpsample,[floor((im.nFFTy-size(sarData,1))/2) 0],0,'pre');
sarDataUpsample = padarray(sarDataUpsample,[ceil((im.nFFTy-size(sarData,1))/2) 0],0,'post');
sarDataUpsample = padarray(sarDataUpsample,[0 floor((im.nFFTx-size(sarData,2))/2)],0,'pre');
sarDataUpsample = padarray(sarDataUpsample,[0 ceil((im.nFFTx-size(sarData,2))/2)],0,'post');


sarDataFFT = fftshift(fftshift(fft(fft((sarDataUpsample),im.nFFTy,1),im.nFFTx,2),1),2)/im.nFFTx/im.nFFTy;

kSx = 2*pi/sar.xStep_m;
kSy = 2*pi/sar.yStep_m;

kX = linspace(-kSx/2,kSx/2,im.nFFTx);
kY = reshape(linspace(-kSy/2,kSy/2,im.nFFTy),[],1);
kZ = single(sqrt(4 * k.^2 - kX.^2 - kY.^2));

focusingFilter = kZ .* exp(-1j * kZ * abs(ant.tx.z0_m - target.zOffset_m));
focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;

sarImage = single(abs(sum(ifft(ifft(sarDataFFT .* focusingFilter,[],1),[],2),3)))';

% sarImage = flip(sarImage,2);

x_m = sar.xStep_m * (-(im.nFFTx-1)/2 : (im.nFFTx-1)/2);
y_m = sar.yStep_m * (-(im.nFFTy-1)/2 : (im.nFFTy-1)/2);

% Resize Image
if max(im.x_m) > max(x_m)
    warning("WARNING: im.nFFTx is too small to see the image!")
end
if max(im.y_m) > max(y_m)
    warning("WARNING: im.nFFTy is too small to see the image!")
end

% indX = x_m >= (im.x_m(1)) & x_m <= (im.x_m(end));
% indY = y_m >= (im.y_m(1)) & y_m <= (im.y_m(end));

im.sarImage = interpn(x_m,y_m,sarImage,im.x_m,im.y_m,'spline',0);
% im.sarImage = imresize(sarImage(indX,indY),[im.numX,im.numY]);

figure
mesh(im.x_m,im.y_m,im.sarImage','FaceColor','interp')
xlabel("x (m)")
ylabel("y (m)")
xlim([im.x_m(1),im.x_m(end)])
ylim([im.y_m(1),im.y_m(end)])
title("Reconstructed Image");
view(2)
