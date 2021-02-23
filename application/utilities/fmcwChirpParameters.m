classdef fmcwChirpParameters
    properties
        f0
        K
        IdleTime_s
        TXStartTime_s
        ADCStartTime_s
        ADCSamples
        ADCSampleTime_time
        ADCSampleTime_sample
        fS
        RampEndTime_s
        fC
        B_max
        B_total
        B_useful
        c
        rangeMax_m
        rangeResolution_m
        k
        lambda_m
    end
    methods
        function obj = fmcwChirpParameters(app)
            obj = update(obj,app);
        end
        
        function obj = update(obj,app)
            obj = getChirpParameters(obj,app);
            obj = computeChirpParameters(obj,app);
            displayChirpParameters(obj,app);
        end
        
        function obj = getChirpParameters(obj,app)
            obj.f0 = app.StartingFrequencyGHzEditField.Value*1e9;
            obj.K = app.FreqSlopeMHzusEditField.Value*1e12;
            obj.IdleTime_s = app.IdleTimeusEditField.Value*1e-6;
            obj.TXStartTime_s = app.TXStartTimeusEditField.Value*1e-6;
            obj.ADCStartTime_s = app.ADCStartTimeusEditField.Value*1e-6;
            obj.ADCSamples = app.ADCSamplesEditField.Value;
            obj.fS = app.SamplingFrequencykspsEditField.Value*1e3;
            obj.RampEndTime_s = app.RampEndTimeusEditField.Value*1e-6;
            obj.fC = app.CenterFrequencyGHzEditField.Value*1e9;
            obj.B_max = app.MaximumBandwidthGHzEditField.Value*1e9;
        end
        
        function obj = computeChirpParameters(obj,app)
            obj.ADCSampleTime_time = obj.RampEndTime_s - obj.ADCStartTime_s - obj.TXStartTime_s;
            obj.ADCSampleTime_sample = obj.ADCSamples/obj.fS;
            
            obj.c = physconst('lightspeed');
            
            if obj.ADCSampleTime_sample > obj.ADCSampleTime_time
                uiconfirm(app.UIFigure,"Not enough time to collect " + obj.ADCSamples + " samples at " + obj.fS*1e-3 + " ksps",'Not Enough Time!',...
                    "Options",{'OK'},'Icon','warning');
            end
            obj.B_total = obj.RampEndTime_s*obj.K;
            obj.B_useful = obj.ADCSampleTime_sample*obj.K;
            
            if obj.B_total > obj.B_max || obj.B_useful > obj.B_max
                warning("Bandwidth exceeds maximum!")
                uiconfirm(app.UIFigure,"Bandwidth exceeds maximum! Reduce ADC Samples and/or Ramp End Time",'Not Enough Bandwidth!',...
                    "Options",{'OK'},'Icon','warning');
            end
            
            obj.rangeMax_m = obj.fS*obj.c/(2*obj.K);
            obj.rangeResolution_m = obj.c/(2*obj.B_useful);
            
            f0_temp = obj.f0 + obj.ADCStartTime_s*obj.K; % This is for ADC sampling offset
            f = f0_temp + (0:obj.ADCSamples-1)*obj.K/obj.fS; % wideband frequency
            
            obj.k = 2*pi*f/obj.c;
            obj.lambda_m = obj.c/(obj.fC);
        end
        
        function displayChirpParameters(obj,app)
            app.TotalBandwidthGHzEditField.Value = obj.B_total*1e-9;
            app.UsefulBandwidthGHzEditField.Value = obj.B_useful*1e-9;
            app.WavelengthmmEditField.Value = obj.lambda_m*1e3;
            app.MaximumRangemEditField.Value = obj.rangeMax_m;
            app.RangeResolutionmmEditField.Value = obj.rangeResolution_m*1e3;
        end
        
        function obj = loadChirpParameters(obj,app)
            if ~exist(app.FMCWLoadNameEditField.Value + ".mat",'file')
                uiconfirm(app.UIFigure,"No file called " + app.FMCWLoadNameEditField.Value + ".mat to load",'Cannot Load',...
                    "Options",{'OK'},'Icon','warning');
                warning("Parameters not loaded!");
                return;
            end
            
            loadPathFull = "./saved/fmcwChirpParameters/" + app.FMCWLoadNameEditField.Value + ".mat";
            load(loadPathFull,"savedfmcw");
            
            app.StartingFrequencyGHzEditField.Value = savedfmcw.f0*1e-9;
            app.FreqSlopeMHzusEditField.Value = savedfmcw.K*1e-12;
            app.IdleTimeusEditField.Value = savedfmcw.IdleTime_s*1e6;
            app.TXStartTimeusEditField.Value = savedfmcw.TXStartTime_s*1e6;
            app.ADCStartTimeusEditField.Value = savedfmcw.ADCStartTime_s*1e6;
            app.ADCSamplesEditField.Value = savedfmcw.ADCSamples;
            app.SamplingFrequencykspsEditField.Value = savedfmcw.fS*1e-3;
            app.RampEndTimeusEditField.Value = savedfmcw.RampEndTime_s*1e6;
            app.CenterFrequencyGHzEditField.Value = savedfmcw.fC*1e-9;
            app.MaximumBandwidthGHzEditField.Value = savedfmcw.B_max*1e-9;
            
            obj = update(obj,app);
        end
        
        function saveChirpParameters(obj,app)
            if exist(app.FMCWSaveNameEditField.Value + ".mat",'file')
                selection = uiconfirm(app.UIFigure,'Are you sure you want to overwrite?','Confirm Overwrite',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("Parameters not saved!");
                    return;
                end
            end
            
            savedfmcw = obj;
            savePathFull = "./saved/fmcwChirpParameters/" + app.FMCWSaveNameEditField.Value + ".mat";
            save(savePathFull,"savedfmcw");
        end
    end
end