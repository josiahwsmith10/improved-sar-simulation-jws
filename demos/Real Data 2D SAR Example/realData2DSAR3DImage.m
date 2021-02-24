% Example with Real Data of a Knife Blade
addpath(genpath("../../"))

%% Load in Saved Array and FMCW
% These can be created in the app
%-------------------------------------------------------------------------%
fmcw = load("fmcw_v1").savedfmcw;

%% Create the SAR Scenario
%-------------------------------------------------------------------------%
sar.x_step_m = fmcw.lambda_m/4;
sar.y_step_m = fmcw.lambda_m/4;

%% Set Imaging Parameters
%-------------------------------------------------------------------------%
im.nFFTx = 512;
im.nFFTy = 512;
im.nFFTz = 512;

im.numX = 128;
im.numY = 128;
im.numZ = 128;
im.x_m = linspace(0.3/im.numX-0.15,0.15,im.numX);
im.y_m = linspace(0.3/im.numY-0.15,0.15,im.numY)';
im.z_m = reshape(linspace(0.1+0.4/im.numZ,0.5,im.numZ),1,1,[]);

%% Load Data
%-------------------------------------------------------------------------%
load rectilinearTest3
% sarData is Y x X x k

%% Reconstruct Image
%-------------------------------------------------------------------------%
sar.sarData = permute(sarData,[2,1,3]);
im = uniform_SISO_2D_array_reconstructImage_3D(sar,fmcw,im,true);

%% Show the Image
%-------------------------------------------------------------------------%
im.dBMin = -7;
plotXYZdB(im.pxyz,im.x_m,im.y_m,im.z_m,[],im.dBMin,"Reconstructed Image",12);
view(-30,17)