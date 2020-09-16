% Test Bench for the Non-App Versions of the Utilities
%% Load in Saved Array and FMCW 
% These can be created in the app
%-------------------------------------------------------------------------%
ant = load("AWR1243").savedant;
fmcw = load("fmcw_v1").savedfmcw;
set(0,'DefaultFigureWindowStyle','docked')

ant.tx.z0_m = 0.25;
ant.rx.z0_m = 0.25;
ant = updateant(ant,fmcw);

%% Create the SAR Scenario
%-------------------------------------------------------------------------%
sar.method = "Rectilinear";
sar.numX = 200;
sar.numY = 25;
sar.xStep_m = fmcw.lambda_m/4;
sar.yStep_m = fmcw.lambda_m*2;

sar.thetaMax_deg = 360;
sar.numTheta = 1024;

sar = updatesar(sar,ant);

%% Create target
%-------------------------------------------------------------------------%
target.numTarget = 1;

%% Simulate Echo Signal
%-------------------------------------------------------------------------%
target = updatetarget(target,sar,fmcw);