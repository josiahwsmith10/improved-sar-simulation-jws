%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.

% Copyright 2011 The MathWorks, Inc.


%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
% fv = stlread('femur.stl');
[f,v] = stlread('./STLLibrary/MaleTorso.stl');
v = v(v(:,2)<-10,:)*2;
% plot3(v(:,1),v(:,2),v(:,3),'.')
% v = v + [0 370 0];
v = v + [0 470 0];
% plot3(v(:,1),v(:,2),v(:,3),'.')

ptCloud = pointCloud(v);
% pcshow(ptCloud);

% ptCloudOut = pcdownsample(ptCloud,'gridAverage',2);
ptCloudOut = pcdownsample(ptCloud,'random',0.025);
pcshow(ptCloudOut);
xlabel('x');ylabel('y');zlabel('z');

v = ptCloudOut.Location;
xyzMatrixT = zeros(size(v));
xyzMatrixT(:,1) = v(:,1);
xyzMatrixT(:,2) = v(:,3);
xyzMatrixT(:,3) = v(:,2);
pcshow(xyzMatrixT);
xlabel('x');ylabel('y');zlabel('z');

% clear f ptCloud ptCloudOut v;

%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);