function fig = initializeFigures()
set(0,'DefaultFigureWindowStyle','docked')

% AntAxes
fig.AntAxes.f = figure;
fig.AntAxes.h = handle(axes);

% AntVirtualAxes
fig.AntVirtualAxes.f = figure;
fig.AntVirtualAxes.h = handle(axes);

% SARAxes
fig.SARAxes.f = figure;
fig.SARAxes.h = handle(axes);

% SARVirtualAxes
fig.SARVirtualAxes.f = figure;
fig.SARVirtualAxes.h = handle(axes);

end