function ant = updateant(app)
%% Get local variables
fmcw = app.fmcw;
fig = app.fig;

ant.tx.z0_m = app.TransmitterZmEditField.Value;
ant.rx.z0_m = app.ReceiverZmEditField.Value;

AntAxes = app.AntAxes;
AntVirtualAxes = app.AntVirtualAxes;

%% Get positions from table
ant.tx.xy = table2array(app.TxTable.Data);
ant.rx.xy = table2array(app.RxTable.Data);

%% Get number or Tx, Rx, and Vx
ant.tx.numTx = sum(ant.tx.xy(:,5));
ant.rx.numRx = sum(ant.rx.xy(:,5));
ant.vx.numVx = ant.tx.numTx * ant.rx.numRx;

%% Get only the enabled elements
if ant.tx.numTx > 0
    ant.tx.xy = ant.tx.xy(logical(ant.tx.xy(:,5)),[1,3])*fmcw.lambda_m + ant.tx.xy(logical(ant.tx.xy(:,5)),[2,4])*1e-3;
else
    ant.tx.xy = [];
end

if ant.rx.numRx > 0
    ant.rx.xy = ant.rx.xy(logical(ant.rx.xy(:,5)),[1,3])*fmcw.lambda_m + ant.rx.xy(logical(ant.rx.xy(:,5)),[2,4])*1e-3;
else
    ant.rx.xy = [];
end

ant.vx.xy = [];
ant.vx.dxy = [];

for indTx = 1:ant.tx.numTx
    ant.vx.xy = cat(1,ant.vx.xy,(ant.tx.xy(indTx,:) + ant.rx.xy)/2);
    ant.vx.dxy = cat(1,ant.vx.dxy,ant.tx.xy(indTx,:) - ant.rx.xy);
end

% Setup Tx/Rx xyz spacing for rectilinear SAR
ant.tx.xyz_m = [ant.tx.xy,ant.tx.z0_m*ones(ant.tx.numTx,1)];
ant.rx.xyz_m = [ant.rx.xy,ant.rx.z0_m*ones(ant.rx.numRx,1)];
ant.tx.xyz_m = repmat(ant.tx.xyz_m,ant.rx.numRx,1);
ant.tx.xyz_m = reshape(ant.tx.xyz_m,ant.tx.numTx,ant.rx.numRx,3);
ant.tx.xyz_m = permute(ant.tx.xyz_m,[2,1,3]);
ant.tx.xyz_m = reshape(ant.tx.xyz_m,ant.vx.numVx,3);
ant.rx.xyz_m = repmat(ant.rx.xyz_m,ant.tx.numTx,1);

ant.tx.xyz_m = reshape(ant.tx.xyz_m,ant.vx.numVx,1,3);
ant.rx.xyz_m = reshape(ant.rx.xyz_m,ant.vx.numVx,1,3);
ant.vx.xyz_m = reshape([ant.vx.xy,ant.tx.z0_m*ones(ant.vx.numVx,1)],ant.vx.numVx,1,3);

% Unwrap ant.tx.xyz_m & ant.rx.xyz_m as [numRx,numTx,3]

%% Scatter plot the Tx, Rx, and Vx elements
if ~ant.tx.numTx || ~ant.rx.numRx
    return;
end

h = AntAxes;
hold(h,'off')
scatter(h,ant.tx.xy(:,1)/fmcw.lambda_m,ant.tx.xy(:,2)/fmcw.lambda_m,'xr');
hold(h,'on')
scatter(h,ant.rx.xy(:,1)/fmcw.lambda_m,ant.rx.xy(:,2)/fmcw.lambda_m,'ob');
legend(h,"Tx","Rx")
xlabel(h,"x (\lambda m)")
xlim(h,[min(min(ant.tx.xy(:,1)/fmcw.lambda_m),min(ant.rx.xy(:,1)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,1)/fmcw.lambda_m),max(ant.rx.xy(:,1)/fmcw.lambda_m))+1])
ylim(h,[min(min(ant.tx.xy(:,2)/fmcw.lambda_m),min(ant.rx.xy(:,2)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,2)/fmcw.lambda_m),max(ant.rx.xy(:,2)/fmcw.lambda_m))+1])
ylabel(h,"y (\lambda m)")
title(h,"Physical Array (x-y)")

h = fig.AntAxes.h;
hold(h,'off')
scatter(h,ant.tx.xy(:,1)/fmcw.lambda_m,ant.tx.xy(:,2)/fmcw.lambda_m,'xr');
hold(h,'on')
scatter(h,ant.rx.xy(:,1)/fmcw.lambda_m,ant.rx.xy(:,2)/fmcw.lambda_m,'ob');
legend(h,"Tx","Rx")
xlabel(h,"x (\lambda m)")
xlim(h,[min(min(ant.tx.xy(:,1)/fmcw.lambda_m),min(ant.rx.xy(:,1)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,1)/fmcw.lambda_m),max(ant.rx.xy(:,1)/fmcw.lambda_m))+1])
ylim(h,[min(min(ant.tx.xy(:,2)/fmcw.lambda_m),min(ant.rx.xy(:,2)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,2)/fmcw.lambda_m),max(ant.rx.xy(:,2)/fmcw.lambda_m))+1])
ylabel(h,"y (\lambda m)")
title(h,"Physical Array (x-y)")


h = AntVirtualAxes;
scatter(h,ant.vx.xy(:,1)/fmcw.lambda_m,ant.vx.xy(:,2)/fmcw.lambda_m,'.k');
legend(h,"Vx")
xlabel(h,"x (\lambda m)")
xlim(h,[min(min(ant.tx.xy(:,1)/fmcw.lambda_m),min(ant.rx.xy(:,1)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,1)/fmcw.lambda_m),max(ant.rx.xy(:,1)/fmcw.lambda_m))+1])
ylim(h,[min(min(ant.tx.xy(:,2)/fmcw.lambda_m),min(ant.rx.xy(:,2)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,2)/fmcw.lambda_m),max(ant.rx.xy(:,2)/fmcw.lambda_m))+1])
ylabel(h,"y (\lambda m)")
title(h,"Virtual Array (x-y)")

h = fig.AntVirtualAxes.h;
scatter(h,ant.vx.xy(:,1)/fmcw.lambda_m,ant.vx.xy(:,2)/fmcw.lambda_m,'.k');
legend(h,"Vx")
xlabel(h,"x (\lambda m)")
xlim(h,[min(min(ant.tx.xy(:,1)/fmcw.lambda_m),min(ant.rx.xy(:,1)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,1)/fmcw.lambda_m),max(ant.rx.xy(:,1)/fmcw.lambda_m))+1])
ylim(h,[min(min(ant.tx.xy(:,2)/fmcw.lambda_m),min(ant.rx.xy(:,2)/fmcw.lambda_m))-1,max(max(ant.tx.xy(:,2)/fmcw.lambda_m),max(ant.rx.xy(:,2)/fmcw.lambda_m))+1])
ylabel(h,"y (\lambda m)")
title(h,"Virtual Array (x-y)")

%% Set App variables
app.ant = ant;
app.fig = fig;

%% Update SAR
app.sar = updatesar(app);
end
