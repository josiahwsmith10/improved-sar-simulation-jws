function plotXYdB(imgXY,X,Y,dBMin,xlab,ylab,titlestr,fontSize)
% Organize in meshgrid format
imgYX = imgXY.';

% Normalize Image
imgYX = imgYX/max(imgYX(:));
imgYX_dB = db(imgYX);
clear imgXYZ imgZXY

imgYX_dB(imgYX_dB<dBMin) = dBMin-100;
imgYX_dB(isnan(imgYX_dB)) = dBMin-100;

mesh(X,Y,imgYX_dB,'FaceColor','interp','EdgeColor','none')
view(2)
colormap('jet')
caxis([dBMin 0])

hc = colorbar;
ylabel(hc, 'dB','fontsize',fontSize)

xlabel(xlab,'fontsize',fontSize)
ylabel(ylab,'fontsize',fontSize)
xlim([X(1),X(end)])
ylim([Y(1),Y(end)])
title(titlestr,'fontsize',fontSize)