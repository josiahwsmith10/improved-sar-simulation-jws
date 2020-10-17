% Test Bench for the Non-App Versions of the Utilities
%% Load in Saved Array and FMCW
% These can be created in the app
%-------------------------------------------------------------------------%
addpath(genpath("./"))
ant = load("AWR1243").savedant;
fmcw = load("fmcw_v1").savedfmcw;
fig = initializeFiguresNAU();

%% Create Array 
%-------------------------------------------------------------------------%
ant.tx.z0_m = 0.25;
ant.rx.z0_m = 0.25;
ant = updateantNAU(ant,fmcw,fig);

%% Create the SAR Scenario
%-------------------------------------------------------------------------%
sar.method = "Rectilinear";
sar.numX = 128;
sar.numY = 16;
sar.xStep_m = fmcw.lambda_m/4;
sar.yStep_m = fmcw.lambda_m*2;

sar.thetaMax_deg = 360;
sar.numTheta = 1024;
sar = updatesarNAU(sar,ant,fig);

%% Set Imaging Parameters
%-------------------------------------------------------------------------%
im.nFFTx = 512;
im.nFFTy = 512;

im.numX = 256;
im.numY = 256;
im.x_m = linspace(0.1/im.numX-0.05,0.05,im.numX);
im.y_m = linspace(0.1/im.numY-0.05,0.05,im.numY)';

%% Create target from png
%-------------------------------------------------------------------------%
target.xStep_m = 0.5e-3;
target.yStep_m = 0.5e-3;
target.xOffset_m = -0.0125; %-0.0125,0.0125
target.yOffset_m = 0.025; %-0.025,0.025
target.zOffset_m = 0;
target = gettarget2DpngNAU("cutout1.png",target,im,fig);
showImScenarioNAU(target,sar,fig);

%% Create target from random points
%-------------------------------------------------------------------------%
target.numTargetMax = 64;
target.zOffset_m = 0;
target.o_x = 2e-3;
target.o_y = 2e-3;
target.ampMin = 0.5;
target.ampMax = 0.5;
target = gettarget2DrandNAU(target,im,fig);
showImScenarioNAU(target,sar,fig);

%% Create target from png and random points
%-------------------------------------------------------------------------%
target.numTargetMax = 64;
target.xStep_m = 0.5e-3;
target.yStep_m = 0.5e-3;
target.xOffset_m = -0.0125; %-0.0125,0.0125
target.yOffset_m = 0.025; %-0.025,0.025
target.zOffset_m = 0;
target.o_x = 2e-3;
target.o_y = 2e-3;
target.ampMin = 0.5;
target.ampMax = 0.5;
target = gettarget2DpngANDrandNAU("circle.png",target,im,fig);
showImScenarioNAU(target,sar,fig);

%% Simulate Echo Signal
%-------------------------------------------------------------------------%
target.isAmplitudeFactor = true;
target.isMIMO = false;
sarData = updatetargetNAU(target,sar,fmcw);

%% Reconstruct Image
%-------------------------------------------------------------------------%
% Convert multistatic-to-monostatic
sarData_y_x_k = reshape(permute(sarData,[1,2,4,3,5]),[],sar.numX,fmcw.ADCSamples);
if target.isMIMO
    pc = exp(-1j * reshape(fmcw.k,1,1,[]) .* repmat(ant.tx.xyz_m(:,:,2) - ant.rx.xyz_m(:,:,2),20,1).^2 / (4 * abs(ant.tx.z0_m - target.zOffset_m)));
    sarDataPC = pc .* sarData_y_x_k;
else
    sarDataPC = sarData_y_x_k;
end

% Consider the SISO virtual array
sar.yStep_m = fmcw.lambda_m/4;
sar.numY = 200;

im = uniform_SISO_2D_array_reconstructImage_2DNAU(sarDataPC,target,fmcw,ant,sar,im);