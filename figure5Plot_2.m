clear;
close all;
clc

load('figure5a.mat')

plotimage((1:nt-2)*5.3,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recorda=seis_record;



load('figure5b.mat')

plotimage((1:nt-2)*5.3,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;



figure;plot((1:nt-2)*5.3,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'r','linewidth',2);

hold on;plot((1:nt-2)*5.3,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'g--','linewidth',2)



% hold on;plot(temp-1*10^-5,'b','linewidth',1);

xlabel('time(ms)')
ylabel('Amp')
legend('traditional FD scheme for L-W','proposed FD scheme 1 for L-W')
grid on
axis([ 0 3000 -0.6 1])
