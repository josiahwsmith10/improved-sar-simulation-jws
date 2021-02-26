classdef uniform_Y_SAR_YZ_RMA < handle
    properties
        sarData
        
        nFFTy
        nFFTz
        
        y_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail = false
        isMult2Mono
        
        zRef_m
        k_vec
        z0_m
        yStep_m
        
        fmcw
        ant
        sar
        target
        im
    end
    
    methods
        function obj = uniform_Y_SAR_YZ_RMA(im)
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
            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;
            
            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;
            
            obj.zRef_m = obj.im.zRef_m;
            obj.k_vec = obj.fmcw.k;
            obj.z0_m = obj.ant.z0_m;
            obj.yStep_m = obj.sar.yStep_m;
        end
        
        function verifyParameters(obj)
            obj.isFail = false;
            
            kZU = single(reshape(linspace(0,2*max(obj.k_vec) - 2*max(obj.k_vec)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            y_m_temp = make_x(obj,obj.sar.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
                        
            if max(abs(obj.y_m)) > max(abs(y_m_temp))
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
            if obj.sar.scanMethod ~= "Linear"
                warning("Must use linear SAR scan along the Y-axis to perform Uniform 1-D SAR 2-D RMA image reconstruction method!");
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
                
        function imXYZ_out = computeReconstruction(obj)
            update(obj);
            
            if ~obj.isFail
                if ~obj.ant.isEPC && obj.isMult2Mono
                    mult2mono(obj);
                end
                
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numY,obj.im.numZ));
            end
        end
        
        function reconstruct(obj)
            % sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
            % Zero-Pad Data: s(y,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            clear sarData
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kZU = single(reshape(linspace(0,2*max(k) - 2*max(k)/obj.nFFTz,obj.nFFTz),1,[]));
            dkZU = kZU(2) - kZU(1);
            
            if obj.isGPU
                reset(gpuDevice)
                k = gpuArray(k);
                kY = gpuArray(kY);
                kZU = gpuArray(kZU);
                sarDataPadded = gpuArray(sarDataPadded);
            end
            
            kYU = repmat(kY,[1,obj.nFFTz]);
            kU = single(1/2 * sqrt(kY.^2 + kZU.^2));
            kZ = single(sqrt((4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2)));
            
            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * obj.z0_m);
            if obj.isAmplitudeFactor
                focusingFilter = kZ .* focusingFilter;
            end
            focusingFilter(4 * k.^2 < kY.^2) = 0;
            
            % Compute FFT across Y & X Dimensions: S(kY,k)
            sarDataFFT = fftshift(fft(conj(sarDataPadded),obj.nFFTy,1),1)/obj.nFFTy;
            clear sarDataPadded sarData
            
            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end
            
            % Stolt Interpolation
            try
                sarImageFFT = interpn(kY(:),k(:), sarDataFFT .* focusingFilter ,kYU,kU,'linear',0);
            catch
                sarImageFFT = zeros(size(kU));
                for indkY = 1:size(kU,1)
                    tempS = squeeze(sarDataFFT(indkY,:) .* focusingFilter(indkY,:));
                    kZTemp = squeeze(kZ(indkY,:));
                    [kZTemp_unique,~,ind_c] = uniquetol(kZTemp);
                    tempS = accumarray(ind_c,tempS);
                    if length(kZTemp_unique) > 2
                        sarImageFFT(indkY,:) = gather(interp1(kZTemp_unique,tempS,kZU,'nearest',0));
                    end
                end
            end
            clear sarDataFFT focusingFilter kY k kYU kU kZ kZU            
            
            % Recover Image by IFT: p(y,z)
            sarImage = single(abs(ifftn(sarImageFFT)));
            clear sarImageFFT focusingFilter
            
            % Declare Spatial Vectors
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
            
            [Y,Z] = ndgrid(obj.y_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(y_m_temp(:),z_m_temp(:),sarImage,Y,Z,'linear',0)));
        end
        
        function displayImage(obj)
            displayImage2D(obj.im,obj.im.y_m,obj.im.z_m,"y (m)","z (m)");
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