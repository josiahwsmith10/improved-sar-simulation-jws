fc = 77e9;
fmin = 73e9;
fmax = 80e9;
vp = physconst('lightspeed');
lambda = vp/fc;

%% Design Patch Antenna
patchElement = design(patchMicrostrip,fc);

patchElement.Tilt = 90;
patchElement.TiltAxis = [1 1 0];

myFigure = gcf;
myFigure.Color = 'w';
pattern(patchElement,fc)

%% Create Array
x = 1e-3*(-10:10);
y = 1e-3*(-10:10);
z = 0;
[X,Y,Z] = ndgrid(x,y,z);
ElementPosition = reshape(cat(4,X,Y,Z),[],3).';

H = phased.ConformalArray('Element',patchElement,'ElementPosition',ElementPosition)