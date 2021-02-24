classdef sarTarget
    properties
        isGPU
        isMIMO
        
        isAmplitudeFactor
        
        isLong
        numTargets
        xyz_m
        amp
        R
        
        png
        stl
        rp
        
        sarData
        
        fig
    end
    methods
        function obj = sarTarget(app)            
            obj = verifyGPU(obj,app);
            obj = initializeFigures(obj);
            obj = update(obj,app);
            
            obj.stl.isLoaded = false;
        end
        
        function obj = update(obj,app)
            app.BeatSignalComputedLamp.Color = "red";
            app.ImageReconstructionCompleteLamp.Color = "red";
            obj = getTarget(obj,app);
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if verifyMIMO(obj,app)
                displayTarget(obj,app);
            else
                displayVirtualTarget(obj,app);
            end
        end
        
        function obj = getTarget(obj,app)
            obj.isMIMO = verifyMIMO(obj,app);
            
            obj.isAmplitudeFactor = app.PathLossCheckBox.Value;
            
            obj.xyz_m = [];
            obj.amp = [];
            
            if app.UseTableofTargetsCheckBox.Value
                obj = getTargetTable(obj,app);
            end
            if app.UsePNGFileCheckBox.Value
                obj = getPNGTarget(obj,app);
            end
            if app.UseSTLFileCheckBox.Value
                obj = getSTLTarget(obj,app);
            end
            
            obj.numTargets = size(obj.xyz_m,1);
            obj.xyz_m = single(obj.xyz_m);
            obj.amp = single(obj.xyz_m);
            
            app.TotalNumTargetsEditField.Value = obj.numTargets;
        end
        
        function obj = getTargetTable(obj,app)
            if isempty(app.TargetTable.Data)
                return;
            end
            temp_table = single(table2array(app.TargetTable.Data));
            
            temp_xyz_m = temp_table(:,1:3);
            temp_amp = temp_table(:,4);
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function obj = getPNGTarget(obj,app)
            obj = getPNGParameters(obj,app);
            
            try
                tMat = imread("./saved/pngstl/" + obj.png.fileName);
            catch
                selection = uiconfirm(app.UIFigure,'Would you liked to locate the PNG file?','PNG File Not Found',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("PNG file not loaded!");
                    return;
                end
                
                [filename,pathname] = uigetfile("./saved/pngstl/*.png","Select Desired PNG File");
                
                if filename == 0
                    warning("PNG file not loaded!");
                    return;
                end
                app.PNGFileNameEditField.Value = filename;
                tMat = imread(string(pathname) + string(filename));
            end
            tMat = tMat(:,:,1);
            tMat(tMat<64) = 0;
            tMat(tMat>0) = 1;
            tMat = ~tMat;
            tMat = fliplr(tMat);
            
            [numY,numX] = size(tMat);
            xAxisT = obj.png.xStep_m * (-(numX-1)/2 : (numX-1)/2) + obj.png.xOffset_m;
            yAxisT = obj.png.yStep_m * (-(numY-1)/2 : (numY-1)/2) + obj.png.yOffset_m;
            zAxisT = obj.png.zOffset_m;
            
            [zT,xT,yT] = meshgrid(zAxisT,xAxisT,yAxisT);
            temp_xyz_m = reshape(permute([xT,yT,zT],[1 3 2]),[],3);
            
            indT = rot90(tMat,-1)==true;
            temp_xyz_m = single(temp_xyz_m(indT,:));
            
            temp_pc = pointCloud(temp_xyz_m);
            temp_pc = pcdownsample(temp_pc,'random',1/obj.png.downsampleFactor);
            temp_xyz_m = temp_pc.Location;
            temp_amp = ones(size(temp_xyz_m,1),1)*obj.png.reflectivity;
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
            
            app.NumTargetsEditField.Value = size(temp_xyz_m,1);
        end
        
        function obj = getPNGParameters(obj,app)
            obj.png.fileName = app.PNGFileNameEditField.Value;
            obj.png.xStep_m = app.XPixelSizemEditField.Value;
            obj.png.yStep_m = app.YPixelSizemEditField.Value;
            obj.png.xOffset_m = app.XOffsetmEditField.Value;
            obj.png.yOffset_m = app.YOffsetmEditField.Value;
            obj.png.zOffset_m = app.ZOffsetmEditField.Value;
            obj.png.reflectivity = app.ReflectivityEditField.Value;
            obj.png.downsampleFactor = app.DownsampleFactorEditField.Value;
        end
        
        function obj = getSTLTarget(obj,app)
            if ~obj.stl.isLoaded
                selection = uiconfirm(app.UIFigure,'Would you liked to load the STL file?','STL File Not Yet Loaded',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("STL file not loaded!");
                    return;
                end
                obj = loadSTL(obj,app);
            end
            
            obj = getSTLParameters(obj,app);
            
            temp_xyz_m = single(zeros(size(obj.stl.v)));
            temp_xyz_m(:,1) = obj.stl.v(:,1);
            temp_xyz_m(:,2) = obj.stl.v(:,3);
            temp_xyz_m(:,3) = obj.stl.v(:,2);
            
            temp_xyz_m = temp_xyz_m + [obj.stl.xOffset_m,0,0];
            temp_xyz_m = temp_xyz_m + [0,obj.stl.yOffset_m,0];
            temp_xyz_m = temp_xyz_m + [0,0,obj.stl.zOffset_m];
            temp_xyz_m = temp_xyz_m(temp_xyz_m(:,3)<obj.stl.zCrop_m,:);
            
            temp_pc = pointCloud(temp_xyz_m);
            temp_pc = pcdownsample(temp_pc,'random',1/obj.stl.downsampleFactor);
            temp_xyz_m = temp_pc.Location;
            temp_amp = ones(size(temp_xyz_m,1),1)*obj.stl.reflectivity;
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
            
            app.NumTargetsEditField_2.Value = size(temp_xyz_m,1);
        end
        
        function obj = getSTLParameters(obj,app)
            obj.stl.fileName = app.STLFileNameEditField.Value;
            obj.stl.zCrop_m = app.ZCropmEditField.Value;
            obj.stl.xOffset_m = app.XOffsetmEditField_2.Value;
            obj.stl.yOffset_m = app.YOffsetmEditField_2.Value;
            obj.stl.zOffset_m = app.ZOffsetmEditField_2.Value;
            obj.stl.reflectivity = app.ReflectivityEditField_2.Value;
            obj.stl.downsampleFactor = app.DownsampleFactorEditField_2.Value;
        end
        
        function obj = loadSTL(obj,app)
            app.LoadSTLLamp.Color = "yellow";
            drawnow
            obj = getSTLParameters(obj,app);
            try
                [~,obj.stl.v] = stlread2011("./saved/pngstl/" + obj.stl.fileName);
            catch
                selection = uiconfirm(app.UIFigure,'Would you liked to locate the STL file?','STL File Not Found',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("STL file not loaded!");
                    obj.stl.isLoaded = false;
                    return;
                end
                
                [filename,pathname] = uigetfile("./saved/pngstl/*.stl","Select Desired STL File");
                
                if filename == 0
                    warning("STL file not loaded!");
                    obj.stl.isLoaded = false;
                    return;
                end
                app.STLFileNameEditField.Value = filename;
                [~,obj.stl.v] = stlread2011(string(pathname) + string(filename));
            end
            obj.stl.v = obj.stl.v*1e-3;
            obj.stl.isLoaded = true;
            app.LoadSTLLamp.Color = "green";
        end
        
        function obj = getRandomTarget(obj,app)
            obj = getRandomParameters(obj,app);
            
            temp_x = obj.rp.xMin_m + (obj.rp.xMax_m - obj.rp.xMin_m)*rand(obj.rp.numTargets,1);
            temp_y = obj.rp.yMin_m + (obj.rp.yMax_m - obj.rp.yMin_m)*rand(obj.rp.numTargets,1);
            temp_z = obj.rp.zMin_m + (obj.rp.zMax_m - obj.rp.zMin_m)*rand(obj.rp.numTargets,1);
            
            temp_xyz_m = [temp_x,temp_y,temp_z];
            temp_amp = obj.rp.ampMin + (obj.rp.ampMax - obj.rp.ampMin)*rand(obj.rp.numTargets,1);
            
            app.TargetTable.Data = array2table([temp_xyz_m,temp_amp]);
        end
        
        function obj = getRandomParameters(obj,app)
            obj.rp.numTargets = app.NumTargetsEditField_3.Value;
            obj.rp.xMin_m = app.XMinmEditField.Value;
            obj.rp.xMax_m = app.XMaxmEditField.Value;
            obj.rp.yMin_m = app.YMinmEditField.Value;
            obj.rp.yMax_m = app.YMaxmEditField.Value;
            obj.rp.zMin_m = app.ZMinmEditField.Value;
            obj.rp.zMax_m = app.ZMaxmEditField.Value;
            obj.rp.ampMin = app.ReflectivityMinEditField.Value;
            obj.rp.ampMax = app.ReflectivityMaxEditField.Value;
        end
        
        function obj = computeTarget(obj,app)
            if obj.isGPU
                reset(gpuDevice)
            end
            
            if obj.isMIMO
                % Get distances
                try
                    obj.R = struct;
                    obj.R.tx = pdist2(app.sar.rx.xyz_m,obj.xyz_m);
                    obj.R.rx = pdist2(app.sar.tx.xyz_m,obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = 1;
                    end
                    obj.isLong = false;
                catch
                    obj.isLong = true;
                end
            else
                % Get distances
                try
                    obj.R = pdist2(app.sar.vx.xyz_m,obj.xyz_m);
                    
                    % Get echo signal
                    R_T_plus_R_R = 2*obj.R;
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R).^2;
                    else
                        amplitudeFactor = 1;
                    end
                    obj.isLong = false;
                catch
                    obj.isLong = true;
                end
            end
            
            if obj.isGPU && ~obj.isLong
                amplitudeFactor = gpuArray(amplitudeFactor);
                R_T_plus_R_R = gpuArray(R_T_plus_R_R);
            end
            % Create the progress dialog
            d = uiprogressdlg(app.UIFigure,'Title','Generating Echo Signal',...
                'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
            
            obj.sarData = single(zeros(size(app.sar.tx.xyz_m,1),app.fmcw.ADCSamples));
            
            try
                if obj.isLong
                    error("oops");
                end
                
                % Fast method
                tocs = zeros(1,app.fmcw.ADCSamples);
                for indK = 1:app.fmcw.ADCSamples
                    if d.CancelRequested
                        warning("Beat Signal not Computed!")
                        obj.sarData = 0;
                        return;
                    end
                    tic
                    temp = exp(1j*app.fmcw.k(indK)*R_T_plus_R_R);
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(:,indK) = single(gather(sum(temp,2)));
                    % Update the progress dialog
                    tocs(indK) = toc;
                    d.Value = indK/app.fmcw.ADCSamples;
                    d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,app.fmcw.ADCSamples);
                end
            catch
                d.Title = "Generating Echo Signal Using Slow Method";
                % Always works method
                tocs = single(zeros(1,app.fmcw.ADCSamples*obj.numTargets));
                count = 0;                
                for indSAR = 1:size(app.sar.rx.xyz_m,1)
                    tic
                    if obj.isMIMO
                        obj.R = struct;
                        obj.R.tx = pdist2(app.sar.rx.xyz_m(indSAR,:),obj.xyz_m);
                        obj.R.rx = pdist2(app.sar.tx.xyz_m(indSAR,:),obj.xyz_m);
                        R_T_plus_R_R = obj.R.tx + obj.R.rx;
                        
                        % Amplitude Factor
                        if obj.isAmplitudeFactor
                            amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                        else
                            amplitudeFactor = 1;
                        end
                    else
                        obj.R = pdist2(app.sar.vx.xyz_m(indSAR,:),obj.xyz_m);
                        
                        % Get echo signal
                        R_T_plus_R_R = 2*obj.R;
                        if obj.isAmplitudeFactor
                            amplitudeFactor = obj.amp./(obj.R).^2;
                        else
                            amplitudeFactor = single(1);
                        end
                    end
                    
                    if obj.isGPU
                        amplitudeFactor = gpuArray(amplitudeFactor);
                        R_T_plus_R_R = gpuArray(R_T_plus_R_R);
                    end
                    
                    for indK = 1:app.fmcw.ADCSamples
                        if d.CancelRequested
                            warning("Beat Signal not Computed!")
                            obj.sarData = 0;
                            return;
                        end
                        count = count + 1;
                        temp = exp(1j*app.fmcw.k(indK)*R_T_plus_R_R);
                        if obj.isAmplitudeFactor
                            temp = amplitudeFactor .* temp;
                        end
                        
                        obj.sarData(indSAR,indK) = single(gather(sum(temp,2)));
                        % Update the progress dialog
                        tocs(count) = toc;
                        d.Value = count/(app.fmcw.ADCSamples*size(app.sar.rx.xyz_m,1));
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,count,app.fmcw.ADCSamples*size(app.sar.rx.xyz_m,1));
                    end
                end
            end
            
            delete(d);
            
            % Reshape echo signal
            obj.sarData = reshape(obj.sarData,[app.sar.sarSize,app.fmcw.ADCSamples]);
        end
        
        function outstr = getEstTime(obj,tocs,currentInd,totalInd)
            avgtoc = mean(tocs(1:currentInd))*(totalInd - currentInd);
            hrrem = floor(avgtoc/3600);
            avgtoc = avgtoc - floor(avgtoc/3600)*3600;
            minrem = floor(avgtoc/60);
            avgtoc = avgtoc - floor(avgtoc/60)*60;
            secrem = round(avgtoc);
            
            if hrrem < 10
                hrrem = "0" + hrrem;
            end
            if minrem < 10
                minrem = "0" + minrem;
            end
            if secrem < 10
                secrem = "0" + secrem;
            end
            outstr = hrrem + ":" + minrem + ":" + secrem;
        end
        
        % Plot/figure functions
        function obj = initializeFigures(obj)
            closeFigures(obj);
            
            % AntAxes
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function closeFigures(obj)
            try
                close(obj.fig.f)
            catch
            end
        end        
        
        function displayTarget(obj,app)
            if isempty(obj.xyz_m)
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = app.sar.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = app.sar.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = app.sar.tx.xyz_m(:,1);
            temp2 = app.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = app.sar.tx.xyz_m(:,3);
            temp2 = app.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = app.sar.tx.xyz_m(:,2);
            temp2 = app.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"MIMO Aperture Image Scenario")
            legend(h,"Tx","Rx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function displayVirtualTarget(obj,app)
            h = obj.fig.h;
            hold(h,'off')
            temp = app.sar.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            hold(h,'on')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = app.sar.tx.xyz_m(:,1);
            temp2 = app.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = app.sar.tx.xyz_m(:,3);
            temp2 = app.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = app.sar.tx.xyz_m(:,2);
            temp2 = app.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"Virtual Aperture Image Scenario")
            legend(h,"Vx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        % Save/load functions
        function saveTarget(obj,app)
            if exist(app.TargetSaveNameEditField.Value + ".mat",'file')
                selection = uiconfirm(app.UIFigure,'Are you sure you want to overwrite?','Confirm Overwrite',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("Target not saved!");
                    return;
                end
            end
            
            savedtarget = obj;
            savedtarget.fig = 0;
            savePathFull = "./saved/targets/" + app.TargetSaveNameEditField.Value + ".mat";
            save(savePathFull,"savedtarget");
        end
        
        function obj = loadTarget(obj,app)
            if ~exist(app.TargetLoadNameEditField.Value + ".mat",'file')
                uiconfirm(app.UIFigure,"No file called " + app.TargetLoadNameEditField.Value + ".mat to load",'Cannot Load',...
                    "Options",{'OK'},'Icon','warning');
                warning("Target not loaded!");
                return;
            end
            
            loadPathFull = "./saved/targets/" + app.TargetLoadNameEditField.Value + ".mat";
            load(loadPathFull,"savedtarget");
            
            savedtarget.fig = obj.fig;
            obj = savedtarget;
            
            app.BeatSignalComputedLamp.Color = "red";
            app.ImageReconstructionCompleteLamp.Color = "red";
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if verifyMIMO(obj,app)
                displayTarget(obj,app);
            else
                displayVirtualTarget(obj,app);
            end
        end
        
        function obj = verifyGPU(obj,app)
            if app.UseGPUCheckBox.Value
                try
                    reset(gpuDevice);
                catch
                    obj.isGPU = false;
                    app.UseGPUCheckBox.Value = false;
                    uiconfirm(app.UIFigure,"Unable to locate Nvidia GPU",'No GPU Available',...
                        "Options",{'OK'},'Icon','warning');
                    return;
                end
                obj.isGPU = true;
            else
                obj.isGPU = false;
            end
        end
        
        function tf = verifyMIMO(obj,app)
            if app.MIMOSwitch.Value == "Use MIMO Array"
                tf = true;
            elseif app.MIMOSwitch.Value == "Use EPC Virtual Elements"
                tf = false;
            end
        end
    end
    
    methods(Static)
    end
end