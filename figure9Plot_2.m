clear;
close all;
clc

load('figure9a.mat')
plotimage((1:nt-2)*3.8,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;


load('figure9b.mat')
plotimage((1:nt-2)*3.8,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


figure;plot((1:nt-2)*3.8,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'r','linewidth',2);

hold on;plot((1:nt-2)*3.8,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'g--','linewidth',2)

% hold on;plot((1:nt-2)*3.8,100*(seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75)))-seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75)))),'m','linewidth',2)

% hold on;plot(temp-1*10^-5,'b','linewidth',1);

xlabel('time(ms)')
ylabel('Amp')
legend('traditional FD scheme for L-W','proposed FD scheme 1 for L-W')
grid on
axis([ 0 5300 -0.62 1])
