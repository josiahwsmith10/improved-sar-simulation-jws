classdef reconstructionAlgorithmTemplate
    properties
        sarData
        
        nFFTx
        nFFTy
        nFFTz
        
        x_m
        y_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail
        
        k_vec
        z0_m
        xStep_m
        yStep_m
        
        fmcw
        ant
        sar
        target
        im
    end
    
    methods
        function obj = reconstructionAlgorithmTemplate(im)
            obj.fmcw = obj.im.fmcw;
            obj.ant = obj.im.ant;
            obj.sar = obj.im.sar;
            obj.target = obj.im.target;
            obj.im = im;
            
            getParameters(obj);
        end
        function obj = update(obj)
            obj = getParameters(obj);
            obj = verifyParameters(obj);
            obj = verifyReconstruction(obj);
        end
        
        function obj = getParameters(obj)
            obj.nFFTx = obj.im.nFFTx;
            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;
            
            obj.x_m = obj.im.x_m;
            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.k_vec = obj.fmcw.k;
            obj.z0_m = obj.ant.z0_m;
            obj.xStep_m = obj.sar.xStep_m;
            obj.yStep_m = obj.sar.yStep_m;
        end
        
        function [obj,imXYZ_out] = computeReconstruction(obj)
            obj = update(obj);
            if ~obj.isFail
                obj = reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            end
        end
        
        function obj = reconstruct(obj)
            
        end
    end
end