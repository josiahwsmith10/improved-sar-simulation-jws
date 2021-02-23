function fig = initializeFiguresNAU()
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

% ImScenarioAxes
fig.ImScenarioAxes.f = figure;
fig.ImScenarioAxes.h = handle(axes);

% ImScenarioVirtualAxes
fig.ImScenarioVirtualAxes.f = figure;
fig.ImScenarioVirtualAxes.h = handle(axes);

% Target2D
fig.Target2D.f = figure;
fig.Target2D.h = handle(axes);

% SAR2D
fig.SAR2D.f = figure;
fig.SAR2D.h = handle(axes);