function fmcw = updatefmcw(app)
%% Get new FMCW parameters from the GUI
fmcw.f0 = app.StartingFrequencyGHzEditField.Value*1e9;
fmcw.K = app.FreqSlopeMHzusEditField.Value*1e12;
fmcw.IdleTime_s = app.IdleTimeusEditField.Value*1e-6;
fmcw.TXStartTime_s = app.TXStartTimeusEditField.Value*1e-6;
fmcw.ADCStartTime_s = app.ADCStartTimeusEditField.Value*1e-6;
fmcw.ADCSamples = app.ADCSamplesEditField.Value;
fmcw.fS = app.SamplingFrequencykspsEditField.Value*1e3;
fmcw.RampEndTime_s = app.RampEndTimeusEditField.Value*1e-6;
fmcw.fC = app.CenterFrequencyGHzEditField.Value*1e9;
fmcw.B_max = app.MaximumBandwidthGHzEditField.Value*1e9;

%% Use the parameters to calculate other necessary parameters
fmcw.ADCSampleTime_time = fmcw.RampEndTime_s - fmcw.ADCStartTime_s - fmcw.TXStartTime_s;
fmcw.ADCSampleTime_sample = fmcw.ADCSamples/fmcw.fS;

fmcw.c = physconst('lightspeed');

if fmcw.ADCSampleTime_sample > fmcw.ADCSampleTime_time
    uiconfirm(app.UIFigure,"Not enough time to collect " + fmcw.ADCSamples + " samples at " + fmcw.fS*1e-3 + " ksps",'Not Enough Time!',...
                    "Options",{'OK'},'Icon','warning');
end
fmcw.B_total = fmcw.RampEndTime_s*fmcw.K;
fmcw.B_useful = fmcw.ADCSampleTime_sample*fmcw.K;

if fmcw.B_total > fmcw.B_max || fmcw.B_useful > fmcw.B_max
    warning("Bandwidth exceeds maximum!")
    uiconfirm(app.UIFigure,"Bandwidth exceeds maximum! Reduce ADC Samples and/or Ramp End Time",'Not Enough Bandwidth!',...
                    "Options",{'OK'},'Icon','warning');
end

fmcw.rangeMax_m = fmcw.fS*fmcw.c/(2*fmcw.K);
fmcw.rangeResolution_m = fmcw.c/(2*fmcw.B_useful);

f0 = fmcw.f0 + fmcw.ADCStartTime_s*fmcw.K; % This is for ADC sampling offset
f = f0 + (0:fmcw.ADCSamples-1)*fmcw.K/fmcw.fS; % wideband frequency

fmcw.k = 2*pi*f/fmcw.c;
fmcw.lambda_m = fmcw.c/(fmcw.fC);

%% Display the parameters
app.TotalBandwidthGHzEditField.Value = fmcw.B_total*1e-9;
app.UsefulBandwidthGHzEditField.Value = fmcw.B_useful*1e-9;
app.WavelengthmEditField.Value = fmcw.lambda_m;
app.MaximumRangemEditField.Value = fmcw.rangeMax_m;
app.RangeResolutionmEditField.Value = fmcw.rangeResolution_m;

end