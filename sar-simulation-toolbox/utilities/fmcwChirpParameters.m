classdef fmcwChirpParameters < handle
    properties (SetObservable)
        f0 = 77e9
        K = 100.036e12
        IdleTime_s = 0
        TXStartTime_s = 0
        ADCStartTime_s = 0
        ADCSamples = 79
        ADCSampleTime_time
        ADCSampleTime_sample
        fS = 2000e3
        RampEndTime_s = 39.98e-6
        fC = 79e9
        B_total = 0
        B_useful = 0
        c = 3e8
        rangeMax_m = 0
        rangeResolution_m = 0
        k
        lambda_m = 0
    end
    methods
        function obj = fmcwChirpParameters()
            attachListener(obj);
        end
        
        function getChirpParameters(obj,in)
            obj.f0 = in.f0;
            obj.K = in.K;
            obj.IdleTime_s = in.IdleTime_s;
            obj.TXStartTime_s = in.TXStartTime_s;
            obj.ADCStartTime_s = in.ADCStartTime_s;
            obj.ADCSamples = in.ADCSamples;
            obj.fS = in.fS;
            obj.RampEndTime_s = in.RampEndTime_s;
            obj.fC = in.fC;
        end
        
        function computeChirpParameters(obj)
            obj.ADCSampleTime_time = obj.RampEndTime_s - obj.ADCStartTime_s - obj.TXStartTime_s;
            obj.ADCSampleTime_sample = obj.ADCSamples/obj.fS;
            
            obj.c = physconst('lightspeed');
            
            if obj.ADCSampleTime_sample > obj.ADCSampleTime_time
                warning("Not enough time to collect " + obj.ADCSamples + " samples at " + obj.fS*1e-3 + " ksps");
            end
            obj.B_total = obj.RampEndTime_s*obj.K;
            obj.B_useful = obj.ADCSampleTime_sample*obj.K;
            
            obj.rangeMax_m = obj.fS*obj.c/(2*obj.K);
            obj.rangeResolution_m = obj.c/(2*obj.B_useful);
            
            f0_temp = obj.f0 + obj.ADCStartTime_s*obj.K;        % This is for ADC sampling offset
            f = f0_temp + (0:obj.ADCSamples-1)*obj.K/obj.fS;    % wideband frequency
            
            obj.k = 2*pi*f/obj.c;
            obj.lambda_m = obj.c/(obj.fC);
        end
        
        function loadChirpParameters(obj,loadName)
            if ~exist(loadName + ".mat",'file')
                warning("No file called " + loadName + ".mat to load. Parameters not loaded!");
                return;
            end
            
            loadPathFull = "./saved/fmcwChirpParameters/" + loadName + ".mat";
            load(loadPathFull,"savedfmcw");
            
            getChirpParameters(obj,savedfmcw);
            computeChirpParameters(obj);
        end
        
        function saveChirpParameters(obj,saveName)
            savePathFull = "./saved/fmcwChirpParameters/" + saveName + ".mat";
            if exist(savePathFull,'file')
                str = input('Are you sure you want to overwrite? Y/N: ','s');
                if str ~= 'Y'
                    warning("Parameters not saved!");
                    return;
                end
            end
            
            savedfmcw = obj;
            save(savePathFull,"savedfmcw");
            disp("FMCW chirp parameters saved to: " + savePathFull);
        end
        
        function attachListener(obj)
            addlistener(obj,{'f0','K','IdleTime_s','TXStartTime_s','ADCStartTime_s','ADCSamples','fS','RampEndTime_s','fC'},'PostSet',@fmcwChirpParameters.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            computeChirpParameters(eventData.AffectedObject);
        end
    end
end