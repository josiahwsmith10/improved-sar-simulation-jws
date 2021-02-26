classdef antennaArray < handle
    properties(SetObservable)
        tx = struct('xy_m',[],'xyz_m',[])
        rx = struct('xy_m',[],'xyz_m',[])
        vx = struct('xy_m','xyz_m')
        z0_m = 0
        fmcw
        
        tableTx = [
            0   0   1.5 5   1
            0.5 0   2.5 5   0
            0   0   3.5 5   1]
        tableRx = [
            0   0   0   0   1
            0   0   0.5 0   1
            0   0   1   0   1
            0   0   1.5 0   1]
        
        fig = struct('f',[],'h',[])
        
        isEPC = false
    end
    methods
        function obj = antennaArray(fmcw)
            attachListener(obj);
            obj.fmcw = fmcw;
        end
        
        function obj = computeAntennaArray(obj)
            obj.tx.xy_m = single(obj.tableTx);
            obj.rx.xy_m = single(obj.tableRx);
            
            % Get number or Tx, Rx, and Vx
            obj.tx.numTx = sum(obj.tx.xy_m(:,5));
            obj.rx.numRx = sum(obj.rx.xy_m(:,5));
            obj.vx.numVx = obj.tx.numTx * obj.rx.numRx;
            
            % Get only the enabled elements
            if obj.tx.numTx > 0
                obj.tx.xy_m = obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[1,3])*obj.fmcw.lambda_m + obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[2,4])*1e-3;
            else
                obj.tx.xy_m = [];
            end
            
            if obj.rx.numRx > 0
                obj.rx.xy_m = obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[1,3])*obj.fmcw.lambda_m + obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[2,4])*1e-3;
            else
                obj.rx.xy_m = [];
            end
            
            obj.vx.xy_m = [];
            obj.vx.dxy = [];
            
            for indTx = 1:obj.tx.numTx
                obj.vx.xy_m = cat(1,obj.vx.xy_m,(obj.tx.xy_m(indTx,:) + obj.rx.xy_m)/2);
                obj.vx.dxy = cat(1,obj.vx.dxy,obj.tx.xy_m(indTx,:) - obj.rx.xy_m);
            end
            
            % Setup Tx/Rx xyz spacing
            obj.tx.xyz_m = [obj.tx.xy_m,obj.z0_m*ones(obj.tx.numTx,1)];
            obj.rx.xyz_m = [obj.rx.xy_m,obj.z0_m*ones(obj.rx.numRx,1)];
            obj.tx.xyz_m = repmat(obj.tx.xyz_m,obj.rx.numRx,1);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.tx.numTx,obj.rx.numRx,3);
            obj.tx.xyz_m = permute(obj.tx.xyz_m,[2,1,3]);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,3);
            obj.rx.xyz_m = repmat(obj.rx.xyz_m,obj.tx.numTx,1);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,1,3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,obj.vx.numVx,1,3);
            obj.vx.xyz_m = reshape([obj.vx.xy_m,obj.z0_m*ones(obj.vx.numVx,1)],obj.vx.numVx,1,3);
            
            temp = mean(obj.vx.xyz_m,1);
            temp(3) = 0;
            
            obj.tx.xyz_m = obj.tx.xyz_m - temp;
            obj.rx.xyz_m = obj.rx.xyz_m - temp;
            obj.vx.xyz_m = obj.vx.xyz_m - temp;
        end
        
        % Plot/figure functions
        function initializeFigures(obj)
            closeFigures(obj)
            set(0,'DefaultFigureWindowStyle','docked')
            
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function closeFigures(obj)
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayAntennaArray(obj)
            if isempty(obj.fig.f)
                initializeFigures(obj);
            end
            
            if obj.isEPC
                displayAntennaArrayEPC(obj);
            else
                displayAntennaArrayMIMO(obj);
            end
        end
        
        function displayAntennaArrayMIMO(obj)
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            scatter(h,obj.tx.xy_m(:,1)/obj.fmcw.lambda_m,obj.tx.xy_m(:,2)/obj.fmcw.lambda_m,'xr');
            hold(h,'on')
            scatter(h,obj.rx.xy_m(:,1)/obj.fmcw.lambda_m,obj.rx.xy_m(:,2)/obj.fmcw.lambda_m,'ob');
            legend(h,"Tx","Rx")
            xlabel(h,"x (\lambda m)")
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))+1])
            ylabel(h,"y (\lambda m)")
            title(h,"Physical Array (x-y)")
        end
        
        function displayAntennaArrayEPC(obj)
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            scatter(h,obj.vx.xy_m(:,1)/obj.fmcw.lambda_m,obj.vx.xy_m(:,2)/obj.fmcw.lambda_m,'.k');
            hold(h,'on')
            legend(h,"Vx")
            xlabel(h,"x (\lambda m)")
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))+1])
            ylabel(h,"y (\lambda m)")
            title(h,"Virtual Array (x-y)")
        end
        
        
        % Save/load functions
        function saveAntennaArray(obj,saveName)
            savePathFull = "./saved/antennaArrays/" + saveName + ".mat";
            if exist(savePathFull,'file')
                str = input('Are you sure you want to overwrite? Y/N: ','s');
                if str ~= 'Y'
                    warning("Antenna array not saved!");
                    return;
                end
            end
            
            savedant = obj;
            save(savePathFull,"savedant");
            disp("Anteanna array saved to: " + savePathFull);
        end
        
        function loadAntennaArray(obj,loadName)
            if ~exist(loadName + ".mat",'file')
                warning("No file called " + loadName + ".mat to load. Antenna array not loaded!");
                return;
            end
            
            loadPathFull = "./saved/antennaArrays/" + loadName + ".mat";
            load(loadPathFull,"savedant");
            
            obj.isEPC = savedant.isEPC;
            obj.z0_m = savedant.z0_m;
            obj.tableTx = savedant.tableTx;
            obj.tableRx = savedant.tableRx;
        end
        
        function attachListener(obj)
            addlistener(obj,{'isEPC','z0_m','tableTx','tableRx','fmcw'},'PostSet',@antennaArray.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            computeAntennaArray(eventData.AffectedObject);
        end
    end
end