classdef sarImage < handle
    properties(SetObservable)
        nFFTx = 512
        nFFTy = 512
        nFFTz = 512
        
        numX = 128
        numY = 128
        numZ = 128
        
        x_m
        y_m
        z_m
        
        xMin_m = -0.2
        xMax_m = 0.2
        
        yMin_m = -0.2
        yMax_m = 0.2
        
        zMin_m = 0
        zMax_m = 0.5
        
        method = "-"
        
        imXYZ
        
        fig = struct('f',[],'h',[])
        dBMin = -25
        fontSize = 12
        vSliceIndex
        
        reconstructor = struct("isFail",false)
        isGPU
        isMult2Mono = false
        zRef_m = 0.25
        zSlice_m
        thetaUpsampleFactor = 1
        
        fmcw
        ant
        sar
        target
    end
    methods
        function obj = sarImage(fmcw,ant,sar,target)
            attachListener(obj);
            obj.fmcw = fmcw;
            obj.ant = ant;
            obj.sar = sar;
            obj.target = target;
        end
        
        function update(obj)
            getImagingParameters(obj);
            if obj.method == "-"
                obj.reconstructor.isFail = true;
            end
        end
        
        function computeImage(obj)
            if isempty(obj.target.sarData)
                warning("Must compute beat signal before image reconstruction!")
                return
            end
            
            update(obj);
            
            if obj.method ~= "-" && ~obj.reconstructor.isFail
                    disp("Attempting image reconstruction using " + obj.method + " method.")
                obj.imXYZ = obj.reconstructor.computeReconstruction();
                if ~obj.reconstructor.isFail
                    disp("Done reconstructing image using " + obj.method + " method.")
                    displayImage(obj);
                end
            end
        end
        
        function getImagingParameters(obj)
            generateAxes(obj);
            
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    obj.reconstructor = uniform_XY_SAR_XYZ_RMA(obj);
                    
                case "Uniform 1-D SAR 2-D RMA"
                    obj.reconstructor = uniform_Y_SAR_YZ_RMA(obj);
                    
                case "Uniform 2-D SAR 2-D FFT"
                    obj.reconstructor = uniform_XY_SAR_XY_FFT(obj);
                    
                case "Uniform 1-D SAR 1-D FFT"
                    obj.reconstructor = uniform_Y_SAR_Y_FFT(obj);
                    
                case "2-D SAR 3-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XYZ_BPA(obj);
                    
                case "2-D SAR 2-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XY_BPA(obj);
                    
                case "Uniform 2-D CSAR 3-D PFA"
                    obj.reconstructor = uniform_thetaY_CSAR_XYZ_PFA(obj);
                    
                case "2-D CSAR 3-D BPA"
                    obj.reconstructor = nonuniform_thetaY_CSAR_XYZ_BPA(obj);
                    
                case "Uniform 1-D CSAR 2-D PFA"
                    obj.reconstructor = uniform_theta_CSAR_XZ_PFA(obj);
                    
                case "1-D CSAR 2-D BPA"
                    obj.reconstructor = nonuniform_theta_CSAR_XZ_BPA(obj);
            end
        end
        
        function generateAxes(obj)
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
        end
        
        % Plot/figure functions
        function initializeFigures(obj)
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
        
        function displayImage(obj)
            if isempty(obj.fig.f)
                initializeFigures(obj);
            end
            
            clf(obj.fig.h)
            obj.reconstructor.displayImage();
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
        
        function attachListener(obj)
            addlistener(obj,{'method','fmcw','ant','sar','target'},'PostSet',@sarImage.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            update(eventData.AffectedObject);
        end
    end
end