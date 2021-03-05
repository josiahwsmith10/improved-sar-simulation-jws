classdef reconstructionAlgorithmTemplate < handle
    % reconstructionAlgorithmTemplate is a reconstructor class that
    % performs the reconstruction algorithm written by the user. This class
    % is meant to be used for easily prototyping new reconstruction
    % algorithms using the FMCW MIMO-SAR Image Reconstruction Toolbox API
    
    properties
        sarData             % Computed beat signal
        
        nFFTx = 512         % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512         % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms
        
        x_m                 % Reconstructed image x axis
        y_m                 % Reconstructed image y axis
        z_m                 % Reconstructed image z axis
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        
        k_vec               % Instantaneous wavenumber vector
        z0_m                % Location of the antenna array in the z-plane
        xStep_m = 1e-3      % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters
        
        fmcw                % fmcwChirpParameters object
        ant                 % sarAntennaArray object
        sar                 % sarScenario object
        target              % sarTarget object
        im                  % sarImage object
    end
    
    methods
        function obj = reconstructionAlgorithmTemplate(im)
            % Set the properties corresponding to the object handles for
            % the imaging scenario and get the parameters from those object
            % handles
            
            obj.fmcw = obj.im.fmcw;
            obj.ant = obj.im.ant;
            obj.sar = obj.im.sar;
            obj.target = obj.im.target;
            obj.im = im;
            
            getParameters(obj);
        end
        
        function update(obj)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the parameters
            
            getParameters(obj);
            verifyParameters(obj);
            verifyReconstruction(obj);
        end
        
        function verifyParameters(obj)
            % Add functionality to verify the parameters of your
            % reconstruction algorithm.
            
            obj.isFail = false;
            
            if (false == true)
                warning("Parameter is wrong!")
                obj.isFail = true;
            end
        end
        
        function verifyReconstruction(obj)
            % Add functionality to verify that the reconstruction can
            % continue (e.g. check that the proper scan has been performed)
            
            obj.isFail = false;
            
            if (false == true)
                warning("Cannot reconstruct!")
                obj.isFail = true;
            end
        end
        
        function getParameters(obj)
            % Get the parameters from the object handles
            
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
        
        function imXYZ_out = computeReconstruction(obj)
            % Update the reconstruction algorithm and attempt the
            % reconstruction
            
            obj = update(obj);
            if ~obj.isFail
                obj = reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            end
        end
        
        function obj = reconstruct(obj)
            % Write your custom image reconstruction algorithm here
            % - Use obj.sarData as the beat signal whose size depends on
            % the scanning method and parameters
        end
        
        function displayImage(obj)
            % Display your reconstructed image using displayImage3D(), displayImage2D(), or displayImage1D()
        end
    end
end