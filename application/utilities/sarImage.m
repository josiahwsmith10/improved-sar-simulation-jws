classdef sarImage
    properties
        nFFTx
        nFFTy
        nFFTz
        
        numX
        numY
        numZ
        
        resizeX
        resizeY
        resizeZ
        
        x_m
        y_m
        z_m
        
        xMin_m
        xMax_m
        
        yMin_m
        yMax_m
        
        zMin_m
        zMax_m
        
        method
        
        imXYZ
        
        fig
        dBMin
        fontSize
        vSliceIndex
        
        reconstructor
        isGPU
    end
    methods
        function obj = sarImage(app)
            obj = initializeFigures(obj);
        end
        
        function obj = update(obj,app)
            obj = getImagingParameters(obj,app);
            if app.ReconstructionAlgorithmDropDown.Value ~= "-"
                obj.reconstructor = obj.reconstructor.update(app,obj);
            else
                obj.reconstructor.isFail = true;
            end
            
            if obj.reconstructor.isFail
                app.ReconstructionAlgorithmDropDown.Value = "-";
            end
        end
        
        function obj = computeImage(obj,app)
            if isempty(app.target.sarData)
                uiconfirm(app.UIFigure,"Must compute beat signal before image reconstruction!",'Beat Signal Error!',...
                    "Options",{'OK'},'Icon','warning');
            end
            
            obj = update(obj,app);
            
            if app.ReconstructionAlgorithmDropDown.Value ~= "-"
                [obj.reconstructor,obj.imXYZ] = obj.reconstructor.computeReconstruction(app,obj);
                
                displayImage(obj,app);
            end            
                
            if isempty(obj.reconstructor.isFail) || obj.reconstructor.isFail
                app.ImageReconstructionCompleteLamp.Color = "red";
            else
                app.ImageReconstructionCompleteLamp.Color = "green";
            end
        end
        
        function obj = getImagingParameters(obj,app)
            obj.method = app.ReconstructionAlgorithmDropDown.Value;
            
            obj.dBMin = app.MindBEditField.Value;
            obj.fontSize = app.FontSizeEditField.Value;
            
            obj.nFFTx = app.NumberofXFFTPointsEditField.Value;
            obj.nFFTy = app.NumberofYFFTPointsEditField.Value;
            obj.nFFTz = app.NumberofZFFTPointsEditField.Value;
            
            obj.xMin_m = app.XMinmEditField_im.Value;
            obj.xMax_m = app.XMaxmEditField_im.Value;
            
            obj.yMin_m = app.YMinmEditField_im.Value;
            obj.yMax_m = app.YMaxmEditField_im.Value;
            
            obj.zMin_m = app.ZMinmEditField_im.Value;
            obj.zMax_m = app.ZMaxmEditField_im.Value;
            
            obj.numX = app.NumXVoxelsEditField.Value;
            obj.numY = app.NumYVoxelsEditField.Value;
            obj.numZ = app.NumZVoxelsEditField.Value;
            
            obj = generateAxes(obj);
            
            obj.isGPU = app.UseGPUCheckBox_2.Value;
            
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    obj.reconstructor = uniform_XY_SAR_XYZ_RMA(app,obj);
                    
                case "Uniform 1-D SAR 2-D RMA"
                    obj.reconstructor = uniform_Y_SAR_YZ_RMA(app,obj);
                    
                case "Uniform 2-D SAR 2-D FFT"
                    obj.reconstructor = uniform_XY_SAR_XY_FFT(app,obj);
                    
                case "Uniform 1-D SAR 1-D FFT"
                    obj.reconstructor = uniform_Y_SAR_Y_FFT(app,obj);
                    
                case "2-D SAR 3-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XYZ_BPA(app,obj);
                    
                case "2-D SAR 2-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XY_BPA(app,obj);
                    
                case "Uniform 2-D CSAR 3-D PFA"
                    obj.reconstructor = uniform_thetaY_CSAR_XYZ_PFA(app,obj);
                    
                case "2-D CSAR 3-D BPA"
                    obj.reconstructor = nonuniform_thetaY_CSAR_XYZ_BPA(app,obj);
                    
                case "Uniform 1-D CSAR 2-D PFA"
                    obj.reconstructor = uniform_theta_CSAR_XZ_PFA(app,obj);
                    
                case "1-D CSAR 2-D BPA"
                    obj.reconstructor = nonuniform_theta_CSAR_XZ_BPA(app,obj);
            end
        end
        
        function obj = generateAxes(obj)
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
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
        
        function displayImage(obj,app)
            obj.dBMin = app.MindBEditField.Value;
            obj.fontSize = app.FontSizeEditField.Value;
            
            clf(obj.fig.h)
            obj.reconstructor.displayImage(obj);
        end
        
        function openInVolumeViewer(obj)
            if ismatrix(obj.imXYZ)
                warning("Cannot open 2D image in volume viewer!");
            end
            
            try
                vvim = abs(permute(obj.imXYZ,[2,1,3]));
                volumeViewer(vvim);
            catch
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
end