%% Make FMCW Parameters
%-------------------------------------------------------------------------%
fmcw.f0 = 77e9;
fmcw.K = 100.036e12;
fmcw.IdleTime_s = 10e-6;
fmcw.TXStartTime_s = 0e-6;
fmcw.ADCStartTime_s = 0e-6;
fmcw.ADCSamples = 79;
fmcw.fS = 2000e3;
fmcw.RampEndTime_s = 39.98e-6;
fmcw.B_max = 4e9;

Params.fmcw = createfmcw(fmcw);
clear fmcw

%% Make Antenna Parameters
%-------------------------------------------------------------------------%
ant.tx.deltaY_m = Params.fmcw.lambda_m*2;
ant.tx.numY = 2;

ant.deltaTxRx = 5e-3;

ant.rx.deltaY_m = Params.fmcw.lambda_m/2;
ant.rx.numY = 4;

sar.numY = 8;
sar.deltaY_m = Params.fmcw.lambda_m*2;

sar = createsar(sar,ant);

function fmcw = createfmcw(fmcw)
fmcw.ADCSampleTime_time = fmcw.RampEndTime_s - fmcw.ADCStartTime_s - fmcw.TXStartTime_s;
fmcw.ADCSampleTime_sample = fmcw.ADCSamples/fmcw.fS;

fmcw.c = physconst('lightspeed');

if fmcw.ADCSampleTime_sample > fmcw.ADCSampleTime_time
    warning("Not enough time to collect " + fmcw.adcSample + " samples at " + fmcw.fS*1e-3 + " ksps",'ERROR');
end
fmcw.B_total = fmcw.RampEndTime_s*fmcw.K;
fmcw.B_useful = fmcw.ADCSampleTime_sample*fmcw.K;

if fmcw.B_total > fmcw.B_max
    warning("Bandwidth exceeds maximum!")
end

fmcw.rangeMax_m = fmcw.fS*fmcw.c/(2*fmcw.K);
fmcw.rangeResolution_m = fmcw.c/(2*fmcw.B_useful);

f0 = fmcw.f0 + fmcw.ADCStartTime_s*fmcw.K; % This is for ADC sampling offset
f = f0 + (0:fmcw.ADCSamples-1)*fmcw.K/fmcw.fS; % wideband frequency

fmcw.k = 2*pi*f/fmcw.c;
fmcw.lambda_m = fmcw.c/(fmcw.f0 + fmcw.B_max/2);
end

function sar = createsar(sar,ant)
ant.tx.locY_m = (0:ant.tx.numY-1)*ant.tx.deltaY_m;
ant.rx.locY_m = ant.tx.locY_m(end) + ant.deltaTxRx + (0:ant.rx.numY-1)*ant.rx.deltaY_m;

ant.vx.locY_m = [];
for indTx = 1:ant.tx.numY
    ant.vx.locY_m = [ant.vx.locY_m,(ant.tx.locY_m(indTx) + ant.rx.locY_m)/2];
end

sar.tx.locY_m = [];
sar.rx.locY_m = [];
sar.vx.locY_m = [];

for indY = 1:sar.numY
    sar.tx.locY_m = [sar.tx.locY_m,ant.tx.locY_m + sar.deltaY_m*(indY-1)];
    sar.rx.locY_m = [sar.rx.locY_m,ant.rx.locY_m + sar.deltaY_m*(indY-1)];
    sar.vx.locY_m = [sar.vx.locY_m,ant.vx.locY_m + sar.deltaY_m*(indY-1)];
end

yAvg = (sar.rx.locY_m(end) + sar.tx.locY_m(end))/4;
sar.tx.locY_m = sar.tx.locY_m - yAvg;
sar.rx.locY_m = sar.rx.locY_m - yAvg;

% figure
% scatter(zeros(size(ant.tx.locY_m)),ant.tx.locY_m)
% hold on
% scatter(zeros(size(ant.rx.locY_m)),ant.rx.locY_m);
% 
% figure
% scatter(zeros(size(sar.tx.locY_m)),sar.tx.locY_m)
% hold on
% scatter(zeros(size(sar.rx.locY_m)),sar.rx.locY_m);

end