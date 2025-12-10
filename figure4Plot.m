clear;
close all;
clc

% clip level 0.5


load('figure4a.mat')
seis_recorda=seis_record;
figure;

hold on;plot((1:nt-2)*3.8,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'r','linewidth',2)

load('figure4b.mat')


seis_recordb=seis_record;


hold on;plot((1:nt-2)*3.8,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'k--','linewidth',2)

load('figure4c.mat')


seis_recordc=seis_record;








% 
hold on;plot((1:nt-2)*3.8,seis_recordc((1:nt-2),75)./max(abs(seis_recordc((1:nt-2),75))),'g','linewidth',2)

load('figure4d.mat')

hold on;plot((1:nt-2)*3,seis_record((1:nt-2),75)./max(abs(seis_record((1:nt-2),75))),'b','linewidth',2)


xlabel('time(ms)')
ylabel('Amp')
legend('traditional FD scheme for L-W','proposed FD scheme 1 for L-W','proposed FD scheme 2 for L-W','T-S domain')
grid on
% axis([ 2190 2250 -8*10^-3 6*10^-3])
box on
