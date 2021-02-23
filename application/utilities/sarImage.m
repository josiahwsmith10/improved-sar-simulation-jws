classdef sarImage
    properties
        nFFTx
        nFFTy
        nFFTz
        
        numX
        numY
        numZ
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
    end
    methods
        function obj = sarImage(app)
            obj = initializeFigures(obj);
        end
        
        function obj = update(obj,app)
            obj = getImagingParameters(obj,app);
            obj = verifyParameters(obj);
        end
        
        function obj = computeImage(obj,app)
%             if ~app.sar.isMIMO
%                 temp.UIFigure = app.UIFigure;
%                 temp.fmcw = app.fmcw;
%                 temp.ant = app.ant;
%                 temp.sar = app.sar;
%                 temp.target = app.target;
%                 temp.sar.yStep_m = app.sar.yStep_m/app.ant.vx.numVx;
%                 temp.sar.numY = app.sar.numY*app.ant.vx.numVx;
%             end
%             obj = uniform_SISO_2D_array_reconstructImage_3D(obj,temp);
            
            obj = update(obj,app);
            
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    obj = uniform_SISO_2D_array_reconstructImage_3D(obj,app);
                case "2-D SAR 3-D BPA"
                    
                case "2-D SAR 2-D BPA"
                    
                case "Uniform 1-D SAR 2-D RMA"
                    
                case "1-D SAR 2-D BPA"
                    
                case "1-D SAR 1-D BPA"
            end
            
            
            
            displayImage(obj);
            app.ImageReconstructionCompleteLamp.Color = "green";
        end
        
        function obj = verifyParameters(obj)
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    
                case "2-D SAR 3-D BPA"
                    
                case "2-D SAR 2-D BPA"
                    
                case "Uniform 1-D SAR 2-D RMA"
                    
                case "1-D SAR 2-D BPA"
                    
                case "1-D SAR 1-D BPA"
            end
        end
        
        function obj = getImagingParameters(obj,app)
            obj.method = app.ReconstructionAlgorithmDropDown.Value;
            
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
        end
        
        function obj = generateAxes(obj)
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
        
            
            x_m = linspace(-0.2,0.2-0.4/obj.numX,obj.numX);
            y_m = linspace(-0.2,0.2-0.4/obj.numY,obj.numY);
            z_m = linspace(0.8,1.2-0.4/obj.numZ,obj.numZ);
        end
        
        % Plot/figure functions
        function obj = initializeFigures(obj)
            closeFigures(obj);
            
            % AntAxes
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
            
            % TODO: temp
            obj.numX= 64;
            obj.numY = 64;
            obj.numZ = 64;
            obj.x_m = linspace(-0.2,0.2-0.4/obj.numX,obj.numX);
            obj.y_m = linspace(-0.2,0.2-0.4/obj.numY,obj.numY);
            obj.z_m = linspace(0.8,1.2-0.4/obj.numZ,obj.numZ);
            
            obj.nFFTx = 512;
            obj.nFFTy = 500;
            obj.nFFTz = 1024;
            
            obj.dBMin = -25;
            obj.fontSize = 12;
        end
        
        function closeFigures(obj)
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayImage(obj)
            h = obj.fig.h;
            if isempty(obj.vSliceIndex)
                obj.vSliceIndex = 1:length(obj.z_m);
            end
            
            % Organize in meshgrid format
            imgZXY = permute(obj.imXYZ,[3,1,2]);
            
            U = reshape(obj.x_m,1,[],1);
            V = reshape(obj.z_m,[],1,1);
            W = reshape(obj.y_m,1,1,[]);
            
            [meshu,meshv,meshw] = meshgrid(U,V,W);
            
            % Normalize Image
            imgZXY = imgZXY/max(imgZXY(:));
            imgZXY_dB = db(imgZXY);
            clear imgXYZ imgZXY
            
            imgZXY_dB(imgZXY_dB<obj.dBMin) = -1e10;
            imgZXY_dB(isnan(imgZXY_dB)) = -1e10;
            
            hs = slice(h,meshu,meshv,meshw,imgZXY_dB,single([]),V(obj.vSliceIndex),single([]));
            set(hs,'FaceColor','interp','EdgeColor','none');
            axis(h,'vis3d');
            
            for kk=1:length(obj.vSliceIndex)
                set(hs(kk),'AlphaData',squeeze(imgZXY_dB(kk+obj.vSliceIndex(1)-1,:,:)),'FaceAlpha','interp');
            end
            
            colormap(h,'jet')
            hc = colorbar(h);
            
            view(h,3)
            daspect(h,[1 1 1])
            caxis(h,[obj.dBMin 0])
            
            ylabel(hc, 'dB','fontsize',obj.fontSize)
            xlabel(h,'x (m)','fontsize',obj.fontSize)
            ylabel(h,'z (m)','fontsize',obj.fontSize)
            zlabel(h,'y (m)','fontsize',obj.fontSize)
            xlim(h,[obj.x_m(1),obj.x_m(end)])
            ylim(h,[obj.z_m(1),obj.z_m(end)])
            zlim(h,[obj.y_m(1),obj.y_m(end)])
            title(h,"Reconstructed Image",'fontsize',obj.fontSize)
        end
        
        function openInVolumeViewer(obj)
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