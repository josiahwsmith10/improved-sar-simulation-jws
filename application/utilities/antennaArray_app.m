classdef antennaArray_app < handle
    properties
        tx
        rx
        vx
        fmcw
        
        % For save/load
        TxTable
        RxTable
        
        fig
        
        isMIMO
    end
    methods
        function obj = antennaArray_app(app)
            obj = initializeFigures(obj,app);
            obj = update(obj,app);
        end
        
        function obj = update(obj,app)
            obj = getAntennaArray(obj,app);
            obj = computeAntennaArray(obj);
            if verifyMIMO(obj,app)
                displayAntennaArray(obj);
            else
                displayVirtualAntennaArray(obj);
            end
        end
        
        function obj = getAntennaArray(obj,app)
            obj.isMIMO = verifyMIMO(obj,app);
            obj.fmcw = app.fmcw;
            obj.tx.z0_m = app.TransmitterZmEditField.Value;
            obj.rx.z0_m = app.ReceiverZmEditField.Value;
            obj.tx.xy_m = table2array(app.TxTable.Data);
            obj.rx.xy_m = table2array(app.RxTable.Data);
        end
        
        function obj = computeAntennaArray(obj)
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
            
            % Setup Tx/Rx xyz spacing for rectilinear SAR
            obj.tx.xyz_m = [obj.tx.xy_m,obj.tx.z0_m*ones(obj.tx.numTx,1)];
            obj.rx.xyz_m = [obj.rx.xy_m,obj.rx.z0_m*ones(obj.rx.numRx,1)];
            obj.tx.xyz_m = repmat(obj.tx.xyz_m,obj.rx.numRx,1);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.tx.numTx,obj.rx.numRx,3);
            obj.tx.xyz_m = permute(obj.tx.xyz_m,[2,1,3]);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,3);
            obj.rx.xyz_m = repmat(obj.rx.xyz_m,obj.tx.numTx,1);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,1,3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,obj.vx.numVx,1,3);
            obj.vx.xyz_m = reshape([obj.vx.xy_m,obj.tx.z0_m*ones(obj.vx.numVx,1)],obj.vx.numVx,1,3);
            
            temp = mean(obj.vx.xyz_m,1);
            temp(3) = 0;
            
            obj.tx.xyz_m = obj.tx.xyz_m - temp;
            obj.rx.xyz_m = obj.rx.xyz_m - temp;
            obj.vx.xyz_m = obj.vx.xyz_m - temp;
        end
        
        % Plot/figure functions
        function obj = initializeFigures(obj,app)
            closeFigures(obj)
            
            try
                obj.fig.f = app.sar.fig.f;
                obj.fig.anth = app.sar.fig.anth;
                obj.fig.sarh = app.sar.fig.sarh;
                if ~(ishandle(obj.fig.anth) && ishandle(obj.fig.sarh))
                    error("error");
                end
            catch
                obj.fig.f = figure;
                obj.fig.anth = handle(subplot(121));
                obj.fig.sarh = handle(subplot(122));
            end
        end
        
        function closeFigures(obj)
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayAntennaArray(obj)
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.anth;
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
        
        function displayVirtualAntennaArray(obj)
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.anth;
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
        function saveAntennaArray(obj,app)
            if exist(app.ArraySaveNameEditField.Value + ".mat",'file')
                selection = uiconfirm(app.UIFigure,'Are you sure you want to overwrite?','Confirm Overwrite',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("Antenna array not saved!");
                    return;
                end
            end
            
            obj.TxTable.Data = app.TxTable.Data;
            obj.RxTable.Data = app.RxTable.Data;
            savedant = obj;
            savePathFull = "./saved/antennaArrays/" + app.ArraySaveNameEditField.Value + ".mat";
            save(savePathFull,"savedant");
        end
        
        function obj = loadAntennaArray(obj,app)
            if ~exist(app.ArrayLoadNameEditField.Value + ".mat",'file')
                uiconfirm(app.UIFigure,"No file called " + app.ArrayLoadNameEditField.Value + ".mat to load",'Cannot Load',...
                    "Options",{'OK'},'Icon','warning');
                warning("Antenna array not loaded!");
                return;
            end
            
            loadPathFull = "./saved/antennaArrays/" + app.ArrayLoadNameEditField.Value + ".mat";
            load(loadPathFull,"savedant");
            
            app.TxTable.Data = savedant.TxTable.Data;
            app.RxTable.Data = savedant.RxTable.Data;
            
            obj = update(obj,app);
        end
        
        
        function tf = verifyMIMO(obj,app)
            if app.MIMOSwitch.Value == "Use MIMO Array"
                tf = true;
            elseif app.MIMOSwitch.Value == "Use EPC Virtual Elements"
                tf = false;
            end
        end
    end
end