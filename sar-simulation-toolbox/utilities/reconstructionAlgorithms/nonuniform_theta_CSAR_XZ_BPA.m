classdef nonuniform_theta_CSAR_XZ_BPA < handle
    properties
        sarData
        
        tx_xyz_m
        rx_xyz_m
        vx_xyz_m
        target_xyz_m
        sizeTarget
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail = false
        isMIMO
        
        k_vec
        
        estTime
        
        fmcw
        ant
        sar
        target
        im
    end
    
    methods
        function obj = nonuniform_theta_CSAR_XZ_BPA(im)
            obj.fmcw = im.fmcw;
            obj.ant = im.ant;
            obj.sar = im.sar;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
            obj.estTime.old = inf;
            obj.estTime.count = 0;
        end
        
        function update(obj)
            getParameters(obj);
            verifyReconstruction(obj);
        end
        
        function getParameters(obj)
            obj.tx_xyz_m = obj.sar.tx.xyz_m;
            obj.rx_xyz_m = obj.sar.rx.xyz_m;
            obj.vx_xyz_m = obj.sar.vx.xyz_m;
            
            [X,Y,Z] = ndgrid(obj.im.x_m(:),0,obj.im.z_m(:));
            obj.sizeTarget = [obj.im.numX,obj.im.numZ];
            
            obj.target_xyz_m = [X(:),Y(:),Z(:)];
            
            obj.sarData = reshape(obj.target.sarData,[],obj.fmcw.ADCSamples);
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.k_vec = obj.fmcw.k;
        end
        
        function verifyReconstruction(obj)
            obj.isFail = false;
            
            if obj.sar.scanMethod ~= "Circular"
                warning("Must use 1-D θ Circular CSAR scan to use 1-D CSAR 2-D BPA image reconstruction method!");
                obj.isFail = true;
                return
            end
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
            if obj.isGPU
                reset(gpuDevice);
            end
            
            % Orient vx_xyz x target_xyz x k
            obj.sarData = reshape(obj.sarData,[],1,length(obj.k_vec));
            k = single(reshape(obj.k_vec,1,1,[]));
            
            try
                % Fast way
                fastBPA(obj,k);
            catch
                mediumBPA(obj,k);
            end
            
            try
                obj.imXYZ = reshape(obj.imXYZ,obj.sizeTarget);
            catch
                warning("WHAT!")
            end
        end
        
        function fastBPA(obj,k)
            obj.imXYZ = single(zeros(1,size(obj.target_xyz_m,1)));
            tocs = single(zeros(1,length(k)));
            for indK = 1:length(k)
                tic
                
                if ~obj.ant.isEPC
                    Rt = pdist2(obj.tx_xyz_m,obj.target_xyz_m);
                    Rr = pdist2(obj.rx_xyz_m,obj.target_xyz_m);
                    amplitudeFactor = Rt .* Rr;
                    R_T_plus_R_R = Rt + Rr;
                else
                    R = pdist2(obj.vx_xyz_m,obj.target_xyz_m);
                    R_T_plus_R_R = 2*R;
                    amplitudeFactor = R.^2;
                end
                
                if obj.isGPU
                    R_T_plus_R_R = gpuArray(R_T_plus_R_R);
                end
                
                bpaKernel = gather(exp(-1j*k(indK)*R_T_plus_R_R));
                if obj.isAmplitudeFactor
                    bpaKernel = bpaKernel .* amplitudeFactor;
                end
                obj.imXYZ = obj.imXYZ + sum(obj.sarData(:,:,indK) .* bpaKernel,1);
                % Update the progress dialog
                tocs(indK) = toc;
                disp("Iteration " + indK + "/" + length(k) + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indK,length(k)));
            end
        end
        
        function mediumBPA(obj,k)
            obj.imXYZ = single(zeros(1,size(obj.target_xyz_m,1)));
            tocs = single(zeros(1,2^14));
            count = 0;
            for indTarget = 1:size(obj.target_xyz_m,1)
                for indK = 1:length(k)
                    tic
                    count = count + 1;
                    
                    if ~obj.ant.isEPC
                        Rt = pdist2(obj.tx_xyz_m,obj.target_xyz_m(indTarget,:));
                        Rr = pdist2(obj.rx_xyz_m,obj.target_xyz_m(indTarget,:));
                        amplitudeFactor = Rt .* Rr;
                        R_T_plus_R_R = Rt + Rr;
                    else
                        R = pdist2(obj.vx_xyz_m,obj.target_xyz_m(indTarget,:));
                        R_T_plus_R_R = 2*R;
                        amplitudeFactor = R.^2;
                    end
                    
                    if obj.isGPU
                        R_T_plus_R_R = gpuArray(R_T_plus_R_R);
                    end
                    
                    bpaKernel = gather(exp(-1j*k(indK)*R_T_plus_R_R));
                    if obj.isAmplitudeFactor
                        bpaKernel = bpaKernel .* amplitudeFactor;
                    end
                    obj.imXYZ(indTarget) = obj.imXYZ(indTarget) + sum(obj.sarData(:,:,indK) .* bpaKernel,1);
                    % Update the progress dialog
                    tocs(mod(count,length(tocs))+1) = toc;
                    disp("Iteration " + count + "/" + length(k)*size(obj.target_xyz_m,1) + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indK,length(k)*size(obj.target_xyz_m,1)));
                end
            end
        end
        
        function displayImage(obj)
            displayImage2D(obj.im,obj.im.x_m,obj.im.z_m,"x (m)","z (m)");
        end
        
        function outstr = getEstTime(obj,tocs,currentInd,totalInd)
            avgtoc = mean(tocs(tocs~=0))*(totalInd - currentInd);
            
            if avgtoc > obj.estTime.old && obj.estTime.count < 100
                obj.estTime.count = obj.estTime.count + 1;
                avgtoc = obj.estTime.old;
            else
                obj.estTime.count = 0;
                obj.estTime.old = avgtoc;
            end
            
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
    end
end