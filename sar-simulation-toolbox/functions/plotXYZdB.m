function plotXYZdB(imgXYZ,X,Y,Z,vSliceIndex,dBMin,titlestr,fontSize)
if isempty(vSliceIndex)
    vSliceIndex = 1:length(Z);
end

imgXYZ = double(imgXYZ);
X = double(X);
Y = double(Y);
Z = double(Z);
vSliceIndex = double(vSliceIndex);

f = figure;
h = handle(axes);

% Organize in meshgrid format
imgZXY = permute(imgXYZ,[3,1,2]);

U = reshape(X,1,[],1);
V = reshape(Z,[],1,1);
W = reshape(Y,1,1,[]);

[meshu,meshv,meshw] = meshgrid(U,V,W);

% Normalize Image
imgZXY = imgZXY/max(imgZXY(:));
imgZXY_dB = db(imgZXY);
clear imgXYZ imgZXY

imgZXY_dB(imgZXY_dB<dBMin) = -1e10;
imgZXY_dB(isnan(imgZXY_dB)) = -1e10;

hs = slice(h,meshu,meshv,meshw,imgZXY_dB,[],V(vSliceIndex),[]);
set(hs,'FaceColor','interp','EdgeColor','none');
set(f,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
axis(h,'vis3d');

for kk=1:length(vSliceIndex)
    set(hs(kk),'AlphaData',squeeze(imgZXY_dB(kk+vSliceIndex(1)-1,:,:)),'FaceAlpha','interp');
end

colormap(h,'jet')
hc = colorbar(h);

view(h,3)
daspect(h,[1 1 1])
caxis(h,[dBMin 0])

ylabel(hc, 'dB','fontsize',fontSize)
xlabel(h,'x (m)','fontsize',fontSize)
ylabel(h,'z (m)','fontsize',fontSize)
zlabel(h,'y (m)','fontsize',fontSize)
xlim(h,[X(1),X(end)])
ylim(h,[Z(1),Z(end)])
zlim(h,[Y(1),Y(end)])
title(h,titlestr,'fontsize',fontSize)