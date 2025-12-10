clear;
close all;
clc

% clip level 0.5


load('figure8a.mat')

% plotimage((1:nt-2)*2.8,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')
seis_recorda=seis_record;
dt

load('figure8b.mat')
dt
% load('figure5c3M8.mat')

% plotimage((1:nt-2)*2.8,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')

seis_recordb=seis_record;

% load('figure5c3M8_2.mat')
% plotimage(seis_record)
% seis_recorde=seis_record;


load('figure8c.mat')
dt
% load('figure5c3M8.mat')

% plotimage((1:nt-2)*2.8,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')

seis_recordc=seis_record;

figure;plot((1:nt-2)*2.8,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'g','linewidth',2)

hold on;plot((1:nt-2)*2.8,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'r:','linewidth',2)
hold on;plot((1:nt-2)*2.8,seis_recordc((1:nt-2),75)./max(abs(seis_recordc((1:nt-2),75))),'k-.','linewidth',2)



%  load('figure8d.mat')
%  seis_recordd=seis_recordp;
% 
% hold on;plot((1:nt-2)*2,seis_recordd((1:nt-2),75)./max(abs(seis_recordd((1:nt-2),75))),'m-.','linewidth',2)





% hold on;plot((1:nt-2)*2,seis_recorde((1:nt-2),75)./max(abs(seis_recorde((1:nt-2),75))),'g','linewidth',2)

% hold on;plot(seis_recorde(:,75)./max(abs(seis_recorde(:,75))),'g','linewidth',2)

% figure;
% % plot(seis_recorda(:,75)./max(abs(seis_recorda(:,75))),'b','linewidth',2);
%
% % hold on;
% plot(1000*(seis_recordb(:,75)./max(abs(seis_recordb(:,75)))-seis_recordc(:,75)./max(abs(seis_recordc(:,75)))),'r','linewidth',2)
% % hold on;plot(seis_recordc(:,75)./max(abs(seis_recordc(:,75)))-(temp./max(abs(temp))),'k','linewidth',2)



% hold on;plot(temp-1*10^-5,'b','linewidth',1);

xlabel('time(ms)')
ylabel('Amp')
legend('traditional FD scheme for L-W','proposed FD scheme 1 for L-W','proposed FD scheme 2 for L-W')
grid on
axis([ 0 5300  -0.6 1])
box on

% axis([ 4300 4500  -0.048 0.043])
% axis([ 4900 5050  -0.012 0.02])
% axis([ 3560 3970  -0.28 0.38])