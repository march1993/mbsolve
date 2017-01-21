%function plotSpeedup()
clear all;
close all;

procs = [1 2 4 8 16 32];

execs = [166.163 86.1294 46.4684 24.22 16.413 12.3223 ];

speedup = execs(1) ./ execs;

% plot
figure;
%xlim([procs(1) procs(end)]);
%ylim([procs(1) procs(end)]);
hold on;
plot(procs, speedup);
plot(procs, procs);

%end