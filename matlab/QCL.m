clear all;
close all;

load ../build-openmp/Optica-Basic.mat

t_x = 998;

x = 0:GridPointSize:XDim;
t = 0:TimeStepSize:SimEndTime;

figure;
plot(x, dm11(:, t_x) + dm22(:, t_x) + dm33(:, t_x));
xlim([0, XDim]);
ylabel('Trace as sanity check');

figure;
plot(dm11(1, :) + dm22(1, :) + dm33(1, :));
ylabel('Trace as sanity check');

figure
plot(t, e(1, :));
ylabel('E-Field at left facet');
