classdef uniform_theta_CSAR_XZ_PFA < handle
    properties
        sarData
        
        nFFTx
        nFFTz
        thetaUpsampleFactor = 1
        
        x_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail = false
        
        theta_rad_vec
        k_vec
        R0_m
        xStep_m
        
        fmcw
        ant
        sar
        target
        im
    end
    
    methods
        function obj = uniform_theta_CSAR_XZ_PFA(im)
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
            obj.nFFTz = obj.im.nFFTz;
            
            obj.x_m = obj.im.x_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.thetaUpsampleFactor = obj.im.thetaUpsampleFactor;
            
            obj.theta_rad_vec = obj.sar.theta_rad;
            obj.k_vec = obj.fmcw.k;
            obj.R0_m = obj.ant.z0_m;
            obj.xStep_m = obj.sar.xStep_m;
        end
        
        function verifyParameters(obj)
            obj.isFail = false;
            
            [x_m_temp,z_m_temp] = getTempXZ(obj);
            
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
            
            if max(obj.z_m) > max(z_m_temp)
                warning("zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            if obj.sar.scanMethod ~= "Circular"
                warning("Must use 1-D θ Circular CSAR scan to use 1-D CSAR 2-D PFA image reconstruction method!");
                obj.isFail = true;
                return
            end
            
            % Ensure single element array
            if obj.ant.tx.numTx ~= 1 || obj.ant.rx.numRx ~= 1
                warning("Array must have only 1 Tx and 1 Rx. Please disable necessary elements.");
                obj.isFail = true;
                return
            end
        end
        
        function [x_m_temp,z_m_temp] = getTempXZ(obj)
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            kR = 2*k;
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
            z_m_temp = double(make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz));
        end
        
        function imXYZ_out = computeReconstruction(obj)
            update(obj);
            
            if ~obj.isFail
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numZ));
            end
        end
        
        function reconstruct(obj)          
            % sarData is of size (sar.numTheta, fmcw.ADCSamples)
            % Zero-Pad Data: s(theta,k)
            sarDataPadded = obj.sarData;
            
            theta_rad = single(reshape(obj.theta_rad_vec,[],1));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            
            kR = 2*k;
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            % Declare Uniform Grid for Interpolation
            theta_radU = atan2(kZU,kXU);
            kU = 1/2*sqrt(kXU.^2 + kZU.^2);
            clear kXU kZU
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),[],1));
                [T,K] = ndgrid(theta_radUp(:),reshape(single(1:size(sarDataPadded,2)),1,[]));
                sarDataPadded = interpn(theta_rad(:),single(1:size(sarDataPadded,2))',sarDataPadded,T,K);
                clear T K
                theta_rad = theta_radUp;
            end
            
            % Compute Azimuth Filter H(Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],1));
            clear kR
            
            % Compute Azimuth Filtered Data: p(theta,k) = IFT_Theta[ S(Theta,k) * H*(Theta,k)]
            azimuthFiltered = ifft(fft(sarDataPadded,[],1) .* conj(azimuthFilterFFT),[],1);
            clear azimuthFilterFFT sarDataPadded
            
            if obj.isGPU
                reset(gpuDevice)
                theta_rad = gpuArray(theta_rad);
                k = gpuArray(k);
                azimuthFiltered = gpuArray(azimuthFiltered);
                theta_radU = gpuArray(theta_radU);
                kU = gpuArray(kU);
            end
            
            % Stolt Interpolation
            sarImageFFT = interpn(theta_rad(:),k(:), azimuthFiltered ,theta_radU,kU,'linear',0);
            clear azimuthFiltered kY theta_rad k azimuthfiltered kYU theta_radU kU
            
            if obj.isGPU
                sarImageFFT = gather(sarImageFFT);
                reset(gpuDevice);
                sarImageFFT = gpuArray(sarImageFFT);
            end
            
            % Recover Image by IFT: p(x,z)
            sarImage = single(abs(ifftshift(ifftshift(ifftn(sarImageFFT),1),2)));
            clear sarImageFFT
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            
            [X,Z] = ndgrid(obj.x_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),z_m_temp(:),sarImage,X,Z,'linear',0)));
        end
        
        function displayImage(obj)
            displayImage2D(obj.im,obj.im.x_m,obj.im.z_m,"x (m)","z (m)");
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end        
    end
end