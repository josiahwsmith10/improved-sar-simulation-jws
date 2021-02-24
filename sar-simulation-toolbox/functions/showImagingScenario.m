function showImagingScenario(target,sar)
figure
h = handle(axes);

hold(h,'off')
temp = sar.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
hold(h,'on')
temp = target.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')

xlabel(h,"x (m)")
temp1 = sar.xyz_m(:,1);
temp2 = target.xyz_m(:,1);
xlim(h,[min([min(temp1),min(temp2)])-0.01,max([max(temp1),max(temp2)])+0.01])
ylabel(h,"z (m)")
temp1 = sar.xyz_m(:,3);
temp2 = target.xyz_m(:,3);
ylim(h,[min([min(temp1),min(temp2)])-0.01,max([max(temp1),max(temp2)])+0.25])
zlabel(h,"y (m)")
temp1 = sar.xyz_m(:,2);
temp2 = target.xyz_m(:,2);
zlim(h,[min([min(temp1),min(temp2)])-0.01,max([max(temp1),max(temp2)])+0.01])
title(h,"Imaging Scenario")
legend(h,"Antenna Array","Target")

view(h,3)
daspect(h,[1 1 1])