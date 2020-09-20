function showImScenarioNAU(target,sar,fig)

h = fig.ImScenarioAxes.h;
hold(h,'off')
temp = sar.tx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
hold(h,'on')
temp = sar.rx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
temp = target.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')

xlabel(h,"x (m)")
temp1 = sar.tx.xyz_m(:,1);
temp2 = sar.rx.xyz_m(:,1);
temp3 = target.xyz_m(:,1);
xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
ylabel(h,"z (m)")
temp1 = sar.tx.xyz_m(:,3);
temp2 = sar.rx.xyz_m(:,3);
temp3 = target.xyz_m(:,3);
ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
zlabel(h,"y (m)")
temp1 = sar.tx.xyz_m(:,2);
temp2 = sar.rx.xyz_m(:,2);
temp3 = target.xyz_m(:,2);
zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
title(h,"MIMO Aperture Image Scenario")
legend(h,"Tx","Rx","Target")


h = fig.ImScenarioVirtualAxes.h;
hold(h,'off')
temp = sar.vx.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
hold(h,'on')
temp = target.xyz_m;
scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')

xlabel(h,"x (m)")
temp1 = sar.tx.xyz_m(:,1);
temp2 = sar.rx.xyz_m(:,1);
temp3 = target.xyz_m(:,1);
xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
ylabel(h,"z (m)")
temp1 = sar.tx.xyz_m(:,3);
temp2 = sar.rx.xyz_m(:,3);
temp3 = target.xyz_m(:,3);
ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
zlabel(h,"y (m)")
temp1 = sar.tx.xyz_m(:,2);
temp2 = sar.rx.xyz_m(:,2);
temp3 = target.xyz_m(:,2);
zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
title(h,"Virtual Aperture Image Scenario")
legend(h,"Vx","Target")