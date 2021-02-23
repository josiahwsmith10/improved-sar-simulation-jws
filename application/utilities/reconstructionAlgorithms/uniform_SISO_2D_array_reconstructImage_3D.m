function obj = uniform_SISO_2D_array_reconstructImage_3D(obj,app)
%% Verify Reconstruction Algorithm
[obj,app,failtf] = verifyReconstruction(obj,app);
if failtf
    warning("Reconstruction Failed!")
    return
end

%% Start the Progress Bar
d = uiprogressdlg(app.UIFigure,'Title','Please Wait',...
    'Message','Reconstructing Image using Uniform FFT Method ');

%% sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
sarData = app.target.sarData;
d.Value = 1/10;

%% Zero-Pad Data: s(y,x,k)
sarDataPadded = sarData;
sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(sarData,1))/2) 0],0,'pre');
sarDataPadded = padarray(sarDataPadded,[0 floor((obj.nFFTx-size(sarData,2))/2)],0,'pre');
clear sarData
d.Value = 2/10;

%% Compute Wavenumbers
k = single(reshape(app.fmcw.k,1,1,[]));
L_x = obj.nFFTx * app.sar.xStep_m;
dkX = 2*pi/L_x;
kX = make_kX(dkX,obj.nFFTx);

L_y = obj.nFFTy * app.sar.yStep_m;
dkY = 2*pi/L_y;
kY = make_kX(dkY,obj.nFFTy)';

kZU = single(reshape(linspace(0,2*max(k),obj.nFFTz),1,1,[]));
dkZU = kZU(2) - kZU(1);

if app.target.isGPU
    reset(gpuDevice)
    k = gpuArray(k);
    kX = gpuArray(kX);
    kY = gpuArray(kY);
    kZU = gpuArray(kZU);
    sarDataPadded = gpuArray(sarDataPadded);
end

kYU = repmat(kY,[1,obj.nFFTx,obj.nFFTz]);
kXU = repmat(kX,[obj.nFFTy,1,obj.nFFTz]);
kU = single(1/2 * sqrt(kX.^2 + kY.^2 + kZU.^2));
kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));
d.Value = 3/10;

%% Compute Focusing Filter
focusingFilter = exp(-1j * kZ * app.ant.tx.z0_m);
if app.target.isAmplitudeFactor
    focusingFilter = kZ .* focusingFilter;
end
focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;
d.Value = 4/10;

%% Compute FFT across Y & X Dimensions: S(kY,kX,k)
sarDataFFT = fftshift(fftshift(fft(fft(conj(sarDataPadded),obj.nFFTy,1),obj.nFFTx,2),1),2)/obj.nFFTx/obj.nFFTy;
clear sarDataPadded sarData

if app.target.isGPU
    focusingFilter = gpuArray(focusingFilter);
end
d.Value = 5/10;

%% Stolt Interpolation
try
    sarImageFFT = interpn(kY(:),kX(:),k(:), sarDataFFT .* focusingFilter ,kYU,kXU,kU,'linear',0);
catch
    sarImageFFT = zeros(size(kU));
    for indkY = 1:size(kU,1)
        for indkX = 1:size(kU,2)
            tempS = squeeze(sarDataFFT(indkY,indkX,:) .* focusingFilter(indkY,indkX,:));
            kZTemp = squeeze(kZ(indkY,indkX,:));
            [kZTemp_unique,~,ind_c] = uniquetol(kZTemp);
            tempS = accumarray(ind_c,tempS);
            if length(kZTemp_unique) > 2
                sarImageFFT(indkY,indkX,:) = interp1(kZTemp_unique,tempS,kZU,'nearest',0);
            end
        end
    end
end
clear sarDataFFT focusingFilter kY kX k kYU kXU kU kZ kZU

if app.target.isGPU
    sarImageFFT = gather(sarImageFFT);
    reset(gpuDevice);
    sarImageFFT = gpuArray(sarImageFFT);
end
d.Value = 6/10;

%% Recover Image by IFT: p(y,x,z)
sarImage = single(abs(ifftn(sarImageFFT)));
clear sarDataFFT focusingFilter
d.Value = 7/10;

%% Reorient Image: p(x,y,z)
sarImage = permute(sarImage,[2,1,3]);
d.Value = 8/10;

%% Declare Spatial Vectors
x_m = make_x(app.sar.xStep_m,obj.nFFTx);
y_m = make_x(app.sar.yStep_m,obj.nFFTy);
z_m = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
d.Value = 9/10;

%% Resize Image
if max(obj.x_m) > max(x_m)
    warning("WARNING: im.nFFTx is too small to see the image!")
end
if max(obj.y_m) > max(y_m)
    warning("WARNING: im.nFFTy is too small to see the image!")
end
if max(obj.z_m) > max(z_m)
    warning("WARNING: im.nFFTz is too small to see the image!")
end

[X,Y,Z] = ndgrid(obj.x_m(:),obj.y_m(:),obj.z_m(:));
obj.imXYZ = single(gather(interpn(x_m(:),y_m(:),z_m(:),sarImage,X,Y,Z,'linear',0)));
if app.target.isGPU
    reset(gpuDevice);
end
d.Value = 10/10;
end

%% Functions
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

function [obj,temp_app,failtf] = verifyReconstruction(obj,app)
failtf = true;
temp_app = struct;

% Ensure array is colinear
if max(diff([app.ant.tx.xy_m(:,1);app.ant.rx.xy_m(:,1)])) > 8*eps
    uiconfirm(app.UIFigure,"MIMO array must be colinear. Please disable necessary elements.",'Array Topology Error!',...
        "Options",{'OK'},'Icon','warning');
    return
end

if app.sar.isMIMO 
    % If using MIMO Array
    % Ensure multistatic-to-monostatic approximation is employed
    
else
    % If using EPC Virtual Elements
    % Ensure virtual array is uniform
    if mean(diff(app.ant.vx.xyz_m(:,2),2)) > eps
        uiconfirm(app.UIFigure,"Virtual antenna array is nonuniform! Change antenna positions.",'Array Topology Error!',...
            "Options",{'OK'},'Icon','warning');
        return
    end
    
    % And sar step size is correct
    if app.sar.yStep_m - mean(diff(app.ant.vx.xyz_m(:,2)))*app.ant.vx.numVx > 8*eps
        uiconfirm(app.UIFigure,"SAR step size is incorrect!",'SAR Scenario Error!',...
            "Options",{'OK'},'Icon','warning');
        return
    end
end

% Everything is okay
temp_app.UIFigure = app.UIFigure;
temp_app.fmcw = app.fmcw;
temp_app.ant = app.ant;
temp_app.sar = app.sar;
temp_app.target = app.target;
temp_app.sar.yStep_m = app.sar.yStep_m/app.ant.vx.numVx;
temp_app.sar.numY = app.sar.numY*app.ant.vx.numVx;

failtf = false;
end