clear;
close all;
clc

% clip level 0.5
load('figure3a.mat')


plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recorda=seis_record;


load('figure3b.mat')

plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure3c.mat')

% load('figure5c3M8.mat')

plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;

% load('figure5c3M8_2.mat')
% plotimage(seis_record)
% seis_recorde=seis_record;


load('figure3e.mat')

% load('figure5c3M8.mat')

plotimage((1:nt-2)*2,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recorde=seis_record;


load('figure3d.mat')
% plotimage(    seis_record)
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);




% figure;plot((1:nt-2)*2,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'b','linewidth',2);
% 
% hold on;plot((1:nt-2)*2,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'r','linewidth',2)
% hold on;plot((1:nt-2)*2,seis_recordc((1:nt-2),75)./max(abs(seis_recordc((1:nt-2),75))),'k','linewidth',2)
% 
% hold on;plot((1:nt-2)*2,temp(1:nt-2)./max(abs(temp)),'m','linewidth',2)
% 
% hold on;plot((1:nt-2)*2,seis_recorde((1:nt-2),75)./max(abs(seis_recorde((1:nt-2),75))),'g','linewidth',2)


figure;plot((1:nt-2)*2,seis_recorda((1:nt-2),75),'b','linewidth',2);

hold on;plot((1:nt-2)*2,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*2,seis_recordc((1:nt-2),75),'k--','linewidth',2)

hold on;plot((1:nt-2)*2,temp(1:nt-2),'m-.','linewidth',2)

hold on;plot((1:nt-2)*2,seis_recorde((1:nt-2),75),'g:','linewidth',2)



xlabel('time(ms)')
ylabel('Amp')
legend('FD scheme without L-W','traditional FD scheme for L-W','proposed FD scheme 1 for L-W','Time dispersion elimination','proposed FD scheme 2 for L-W')
grid on
% axis([ 2190 2250 -8*10^-3 6*10^-3])
