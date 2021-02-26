classdef uniform_thetaY_CSAR_XYZ_PFA < handle
    properties
        sarData
        
        nFFTx
        nFFTy
        nFFTz
        thetaUpsampleFactor
        
        x_m
        y_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail = false
        isMult2Mono
        
        zRef_m
        theta_rad_vec
        k_vec
        R0_m
        xStep_m
        yStep_m
        
        fmcw
        ant
        sar
        target
        im
    end
    
    methods
        function obj = uniform_thetaY_CSAR_XYZ_PFA(im)
            obj.fmcw = im.fmcw;
            obj.ant = im.ant;
            obj.sar = im.sar;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
        end
        
        function update(obj)
            getParameters(obj);
            verifyParameters(obj);
            verifyReconstruction(obj);
        end
        
        function getParameters(obj)
            obj.nFFTx = obj.im.nFFTx;
            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;
            obj.thetaUpsampleFactor = obj.im.thetaUpsampleFactor;
            
            obj.x_m = obj.im.x_m;
            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;
            
            obj.zRef_m = obj.im.zRef_m;
            obj.theta_rad_vec = obj.sar.theta_rad;
            obj.k_vec = obj.fmcw.k;
            obj.R0_m = obj.ant.z0_m;
            obj.xStep_m = obj.sar.xStep_m;
            obj.yStep_m = obj.sar.yStep_m;
        end
        
        function verifyParameters(obj)
            obj.isFail = false;
            
            [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj);
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                warning("xMax_m is too large for nFFTx. Decrease xMax_m or increase nFFTx")
                obj.isFail = true;
                return;
            end
            if obj.thetaUpsampleFactor*obj.sar.numTheta > obj.nFFTx
                warning("thetaUpsampleFactor is too large for nFFTx. Decrease thetaUpsampleFactor or increase nFFTx")
                obj.isFail = true;
                return;
            end
            if max(obj.y_m) > max(y_m_temp)
                warning("yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy")
                obj.isFail = true;
                return;
            end
            if max(obj.z_m) > max(z_m_temp)
                warning("zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            if obj.sar.scanMethod ~= "Cylindrical"
                warning("Must use 2-D θY Cylindrical CSAR scan to use 2-D CSAR 3-D PFA image reconstruction method!");
                obj.isFail = true;
                return
            end
            
            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*eps
                warning("MIMO array must be colinear. Please disable necessary elements.");
                obj.isFail = true;
                return
            end
            
            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > eps
                warning("Virtual antenna array is nonuniform! Change antenna positions.");
                obj.isFail = true;
                return
            end
            
            % And sar step size is correct
            if obj.sar.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*eps
                warning("SAR step size is incorrect!");
                obj.isFail = true;
                return
            end
            
            if ~obj.ant.isEPC && ~obj.isMult2Mono
                % If using MIMO Array
                % Ensure multistatic-to-monostatic approximation is employed
                warning("Must use multistatic-to-monostatic approximation to use uniform method!");
                obj.isFail = true;
                return
            end
            
            % Everything is okay to continue
            obj.yStep_m = obj.sar.yStep_m/obj.ant.vx.numVx;
        end
        
        function [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj)
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kR = single(sqrt(4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2));
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),1,[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            x_m_temp = double(make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx));
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            z_m_temp = double(make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz));
        end
        
        function imXYZ_out = computeReconstruction(obj)
            update(obj);
            
            if ~obj.isFail
                if ~obj.ant.isEPC && obj.isMult2Mono
                    mult2mono(obj);
                end
                
%                 reconstruct4segs(obj);
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numY,obj.im.numZ));
            end
        end
        
%         function reconstruct4segs(obj)
%             theta_rad_orig = obj.theta_rad_vec;
%             sarData_orig = obj.sarData;
%             
%             % 1st segment
%             obj.theta_rad_vec = theta_rad_orig(1:end/4);
%             obj.sarData = sarData_orig(:,1:end/4,:);
%             reconstruct(obj);
%             im1 = obj.imXYZ;
%             
%             % 2nd segment
%             obj.theta_rad_vec = theta_rad_orig(end/4+1:end/2);
%             obj.sarData = sarData_orig(:,end/4+1:end/2,:);
%             reconstruct(obj);
%             im2 = obj.imXYZ;
%             
%             % 3rd segment
%             obj.theta_rad_vec = theta_rad_orig(end/2+1:end*3/4);
%             obj.sarData = sarData_orig(:,end/2+1:end*3/4,:);
%             reconstruct(obj);
%             im3 = obj.imXYZ;
%             
%             % 4th segment
%             obj.theta_rad_vec = theta_rad_orig(end*3/4+1:end);
%             obj.sarData = sarData_orig(:,end*3/4+1:end,:);
%             reconstruct(obj);
%             im4 = obj.imXYZ;
%             
%             obj.imXYZ = im1 + im2 + im3 + im4;
%         end
        
        function reconstruct(obj)            
            % sarData is of size (sar.numY, sar.numTheta, fmcw.ADCSamples)
            % Zero-Pad Data: s(y,theta,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kR = single(sqrt(4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2));
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kYU = repmat(kY,[1,obj.nFFTx,obj.nFFTz]);
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),1,[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            % Declare Uniform Grid for Interpolation
            theta_radU = repmat(atan2(kZU,kXU),[obj.nFFTy,1,1]);
            kU = 1/2*sqrt(kXU.^2 + kZU.^2 + kY.^2);
            clear kXU kZU
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),1,[]));
                [Y,T,K] = ndgrid(single(1:size(sarDataPadded,1))',theta_radUp,reshape(single(1:size(sarDataPadded,3)),1,1,[]));
                sarDataPadded = interpn(single(1:size(sarDataPadded,1))',theta_rad(:),single(1:size(sarDataPadded,3))',sarDataPadded,Y,T,K);
                clear Y T K
                theta_rad = theta_radUp;
            end
            
            % Compute Azimuth Filter H(kY,Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],2));
            clear kR
            
            % Compute Azimuth Filtered Data: p(kY,theta,k) = IFT_Theta[ S(kY,Theta,k) * H*(kY,Theta,k)]
            azimuthFiltered = ifft(fftshift(fft(fft(sarDataPadded,[],2),obj.nFFTy,1),1) .* conj(azimuthFilterFFT),[],2);
            clear azimuthFilterFFT sarDataPadded
            
            if obj.isGPU
                reset(gpuDevice)
                kY = gpuArray(kY);
                theta_rad = gpuArray(theta_rad);
                k = gpuArray(k);
                azimuthFiltered = gpuArray(azimuthFiltered);
                kYU = gpuArray(kYU);
                theta_radU = gpuArray(theta_radU);
                kU = gpuArray(kU);
            end
            
            % Stolt Interpolation
            sarImageFFT = interpn(kY(:),theta_rad(:),k(:), azimuthFiltered ,kYU,theta_radU,kU,'linear',0);
            clear azimuthFiltered kY theta_rad k azimuthfiltered kYU theta_radU kU
            
            if obj.isGPU
                sarImageFFT = gather(sarImageFFT);
                reset(gpuDevice);
                sarImageFFT = gpuArray(sarImageFFT);
            end
            
            % Recover Image by IFT: p(y,x,z)
            sarImage = single(ifftshift(ifftshift(ifftn(sarImageFFT),2),3));
            clear sarImageFFT
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            
            [X,Y,Z] = ndgrid(obj.x_m(:),obj.y_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),y_m_temp(:),z_m_temp(:),sarImage,X,Y,Z,'linear',0)));
        end
        
        function displayImage(obj)
            displayImage3D(obj.im);
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end
        
        function kX = make_kX(obj,dkX,nFFTx)
            if mod(nFFTx,2)==0
                kX = dkX * ( -nFFTx/2 : nFFTx/2-1 );
            else
                kX = dkX * ( -(nFFTx-1)/2 : (nFFTx-1)/2 );
            end
            kX = single(kX);
        end
    end
end