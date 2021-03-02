classdef sarScenario < handle
    properties(SetObservable)
        X_m
        Y_m
        Z_m
        xyz_m
        tx
        rx
        vx
        xSize_m
        ySize_m
        xStep_m = 1e-3
        yStep_m = 8e-3
        numX = 1
        numY = 1
        numTheta = 1
        x_m
        y_m
        z_m
        theta_rad
        thetaMax_deg = 360
        
        scanMethod
        sarSize
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
        ant
    end
    methods
        function obj = sarScenario(ant)
            obj.ant = ant;
            attachListener(obj);
        end
        
        function obj = computeSarScenario(obj)
            switch obj.scanMethod
                case "Linear"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyLinear(obj)
                        return
                    end
                    
                    obj = getsarAxes(obj);
                    
                    [obj.X_m,obj.Y_m,obj.Z_m] = ndgrid(obj.x_m,obj.y_m,obj.z_m);
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if ~obj.ant.isEPC
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,3]
                        obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY];
                    else
                        obj.sarSize = obj.ant.vx.numVx*obj.numY;
                    end
                    
                case "Rectilinear"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyRectilinear(obj)
                        return
                    end
                    
                    obj = getsarAxes(obj);
                    
                    [obj.Y_m,obj.X_m,obj.Z_m] = ndgrid(obj.y_m,obj.x_m,obj.z_m);
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if ~obj.ant.isEPC
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numX,3]
                        obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numX];
                    else
                        obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numX];
                    end
                    
                case "Circular"
                    if verifyCircular(obj)
                        return
                    end
                    
                    % Verify single element array
                    if obj.ant.tx.numTx ~= 1 || obj.ant.rx.numRx ~= 1
                        error("Array must have only 1 Tx and 1 Rx. Please disable necessary elements or change the antenna array.")
                    end
                    
                    obj = getsarAxes(obj);
                    
                    obj.x_m = obj.ant.z0_m*cos(obj.theta_rad);
                    obj.y_m = zeros(size(obj.theta_rad));
                    obj.z_m = obj.ant.z0_m*sin(obj.theta_rad) - obj.ant.z0_m;
                    
                    obj.xyz_m = reshape([obj.x_m(:),obj.y_m(:),obj.z_m(:)],1,[],3);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numTheta,3]
                    obj.sarSize = obj.numTheta;
                    
                case "Cylindrical"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyCylindrical(obj)
                        return
                    end
                    
                    obj = getsarAxes(obj);
                    
                    obj.x_m = obj.ant.z0_m*cos(obj.theta_rad);
                    obj.z_m = obj.ant.z0_m*sin(obj.theta_rad) - obj.ant.z0_m;
                    
                    % Use y' as first dimension -> so we obtain s(y',theta)
                    obj.X_m = repmat(obj.x_m(:).',obj.numY,1);
                    obj.Z_m = repmat(obj.z_m(:).',obj.numY,1);
                    obj.Y_m = repmat(obj.y_m(:),1,obj.numTheta);
                    
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if ~obj.ant.isEPC
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numTheta,3]
                        obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numTheta];
                    else
                        obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numTheta];
                    end
            end
            
            obj.tx.xyz_m = single(obj.tx.xyz_m);
            obj.rx.xyz_m = single(obj.rx.xyz_m);
            obj.vx.xyz_m = single(obj.vx.xyz_m);
        end
        
        function obj = getsarAxes(obj)
            % Update number of steps
            obj.xSize_m = obj.numX * obj.xStep_m;
            obj.ySize_m = obj.numY * obj.yStep_m;
            
            % Create synthetic aperture step axes
            obj.x_m = (-(obj.numX - 1)/2 : (obj.numX - 1)/2) * obj.xStep_m;
            obj.y_m = (-(obj.numY - 1)/2 : (obj.numY - 1)/2) * obj.yStep_m;
            obj.theta_rad = linspace(0,obj.thetaMax_deg - obj.thetaMax_deg/obj.numTheta,obj.numTheta)*2*pi/360;
            obj.z_m = 0;
        end
        
        % Figure functions
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
        
        function displaySarScenario(obj)
            if isempty(obj.fig.f)
                initializeFigures(obj);
            end
            
            if obj.ant.isEPC
                displaySarScenarioEPC(obj);
            else
                displaySarScenarioMIMO(obj);
            end
        end
        
        function displaySarScenarioMIMO(obj)
            if obj.scanMethod == "-"
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = obj.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            xlabel(h,"x (m)")
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            title(h,"MIMO Synthetic Aperture")
            legend(h,"Tx","Rx")
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function displaySarScenarioEPC(obj)
            if obj.scanMethod == "-"
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            xlabel(h,"x (m)")
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            title(h,"Virtual Synthetic Aperture")
            legend(h,"Vx")
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function isFail = verifyLinearity(obj)
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*eps
                warning("MIMO array must be colinear. Please disable necessary elements.")
                obj.scanMethod = "-";
                isFail = true;
            else
                isFail = false;
            end
        end
        
        function isFail = verifyLinear(obj)
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for linear scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for linear scan!")
                isFail = true;
            end
            if obj.numTheta ~= 1
                warning("numTheta must be 1 for linear scan!")
                isFail = true;
            end
        end
        
        function isFail = verifyRectilinear(obj)
            isFail = false;
            if obj.numX < 1
                warning("numX must be an integer greater than 1 for rectilinear scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for rectilinear scan!")
                isFail = true;
            end
            if obj.numTheta ~= 1
                warning("numTheta must be 1 for rectilinear scan!")
                isFail = true;
            end            
        end
        
        function isFail = verifyCircular(obj)
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for circular scan!")
                isFail = true;
            end
            if obj.numY ~= 1
                warning("numY must be 1 for circular scan!")
                isFail = true;
            end
            if obj.numTheta < 1
                warning("numTheta must be an integer greater than 1 for circular scan!")
                isFail = true;
            end            
        end
        
        function isFail = verifyCylindrical(obj)
            % Verify proper antenna and SAR scenario parameters for a
            % cylindrical scan
            
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for cylindrical scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for cylindrical scan!")
                isFail = true;
            end
            if obj.numTheta < 1
                warning("numTheta must be an integer greater than 1 for cylindrical scan!")
                isFail = true;
            end            
        end
        
        function attachListener(obj)
            % Attaches the listener to the object handle
            
            addlistener(obj,{'scanMethod','xStep_m','numX','yStep_m','numY','thetaMax_deg','numTheta','ant'},'PostSet',@sarScenario.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            % Recomputes the SAR antenna positions if watched properties are changed
            
            computeSarScenario(eventData.AffectedObject);
        end
    end
end