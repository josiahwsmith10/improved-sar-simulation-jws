c = physconst('lightspeed');

fc = 77e9;
lambda = c/fc;

waveform = phased.FMCWWaveform('SweepTime',32e-9,'SweepBandwidth',4e9,'SweepDirection','Up','SampleRate',4e9);

% ant = phased.IsotropicAntennaElement;

ant = design(patchMicrostrip,fc);
ant.Tilt = 0;
ant.TiltAxis = [0 1 0];

transmitter = phased.Transmitter('PeakPower',1e-3,'Gain',30);

receiver = phased.ReceiverPreamp('Gain',0,'NoiseFigure',4.5,'SampleRate',waveform.SampleRate);

target = phased.RadarTarget('MeanRCS',100,'PropagationSpeed',c,'OperatingFrequency',fc);

channel = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,'SampleRate',waveform.SampleRate,'TwoWayPropagation',true);

tic
% Get echo signal
sig = waveform();
txsig = transmitter(sig);

txsig = channel(txsig,[0;0;0],[0;0;1],[0;0;0],[0;0;0]);
txsig = target(txsig);

txsig = receiver(txsig);
dechirpsig = dechirp(txsig,sig);
toc

Rmax = 4e9 * c / (2*4e9/32e-9);

plot(linspace(0,Rmax,512),abs(fft(dechirpsig,512)))


%% Array
c = physconst('lightspeed');

fc = 77e9;
lambda = c/fc;
T = 32e-9;
B = 4e9;
fs = 4e9;

waveform = phased.FMCWWaveform('SweepTime',T,'SweepBandwidth',B,'SweepDirection','Up','SampleRate',fs);

ant = design(patchMicrostrip,fc);

transmitter = phased.Transmitter('PeakPower',1e-3,'Gain',30);
radiator = phased.Radiator('Sensor',ant,'OperatingFrequency',fc);

receiver = phased.ReceiverPreamp('Gain',30,'NoiseFigure',4.5,'SampleRate',waveform.SampleRate);
collector = phased.Collector('Sensor',ant,'OperatingFrequency',fc);

target.reflector = phased.RadarTarget('MeanRCS',100,'PropagationSpeed',c,'OperatingFrequency',fc);

channel = phased.WidebandFreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);
xp = 1e-3*(-32:31);
yp = 0;
zp = 0;
[Xp,Yp,Zp] = ndgrid(xp,yp,zp);
sar.xyz_m = reshape(cat(4,Xp,Yp,Zp),[],3).';

target.xyz_m = [0;0;01];

% Get echo signal
sig = waveform();
txsig = transmitter(sig);

txsig = radiator(txsig,0);

txsig = channel(repmat(txsig(:),1,size(sar.xyz_m,2)),sar.xyz_m,target.xyz_m,zeros(size(sar.xyz_m)),zeros(size(target.xyz_m)));
txsig = target.reflector(txsig);

txsig = collector(txsig,zeros(2,size(txsig,2)));
txsig = receiver(txsig);
dechirpsig = dechirp(txsig,sig);

Rmax = 4e9 * c / (2*4e9/32e-9);

plot(linspace(0,Rmax,512),abs(fft(dechirpsig,512)))

x = linspace(-0.5,0.5-1/128,128);
z = linspace(0,2-2/128,128);
k = 2*pi/c*(fc + (0:size(dechirpsig,1)-1)*B/T/fs);

im = bpa2D(dechirpsig.',k,x,z,xp,true);

meshyz(im,x,z,"Simualted",2)

%% Functions

function image = bpa2D(syk,k,y,z,yp,isGPU)
image = zeros(length(y),length(z));

yp = reshape(yp,1,1,[]);
k = reshape(k,1,1,1,[]);
syk = reshape(syk,1,1,size(syk,1),size(syk,2));

y = reshape(y,[],1);
z = reshape(z,1,[]);

if isGPU
    syk = gpuArray(single(syk));
    k = gpuArray(single(k));
    yp = gpuArray(single(yp));
    z = gpuArray(single(z));
    image = gpuArray(single(image));
end

try
    for indY = 1:length(y)
        R = sqrt( (yp - y(indY)).^2 + z.^2 );
        bpaKernel = exp(-1j*k*2.*R);
        image(indY,:) = sum(bpaKernel .* syk,[3,4]);
    end
catch
    for indY = 1:length(y)
        for indZ = 1:length(z)
            R = sqrt( (yp - y(indY)).^2 + z(indZ).^2 );
            bpaKernel = exp(-1j*k*2.*R);
            image(indY,indZ) = sum(bpaKernel .* syk,'all');
        end
    end
end
end

function meshyz(im_yz,y,z,titleStr,viewMat,dBmin)
if nargin > 5
    im_yz = db(abs(im_yz)./max(im_yz(:)));
    mesh(y,z,im_yz.','FaceColor','interp','EdgeColor','none');
    xlabel("y (m)")
    ylabel("z (m)")
    title(titleStr)
    view(viewMat)
    zlim([dBmin,0])
    colorbar
    caxis([dBmin,0])
else
    mesh(y,z,abs(im_yz).','FaceColor','interp','EdgeColor','none');
    xlabel("y (m)")
    ylabel("z (m)")
    title(titleStr)
    view(viewMat)
end
end