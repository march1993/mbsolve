clear all;
close all;

load ../build-openmp/Optica-Basic.mat

t_x = 990;

x = 0:GridPointSize:XDim;
%t = 0:TimeStepSize:SimEndTime;

figure;
plot(x, dm11(:, t_x) + dm22(:, t_x) + dm33(:, t_x));
xlim([0, XDim]);
ylabel('Trace as sanity check');