clear;
close all;
clc



load('figure7a.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;


load('figure7b.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;



load('figure7c.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordc=seis_record;



load('figure7d.mat')
% plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')
seis_recordd=seis_record;

NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);








load('figure7e.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorde=seis_record;

figure ;plot((1:nt-2)*1.5,seis_recorda((1:nt-2),75)./max(abs(seis_recorda((1:nt-2),75))),'b','linewidth',2)
hold on;plot((1:nt-2)*1.5,seis_recordb((1:nt-2),75)./max(abs(seis_recordb((1:nt-2),75))),'r','linewidth',2)
hold on;plot((1:nt-2)*1.5,seis_recordc((1:nt-2),75)./max(abs(seis_recordc((1:nt-2),75))),'k-.','linewidth',2)
hold on;plot((1:nt-2)*1.5,temp(1:nt-2)./max(abs(temp)),'m--','linewidth',2)
hold on;plot((1:nt-2)*1.5,seis_recorde((1:nt-2),75)./max(abs(seis_recorde((1:nt-2),75))),'g:','linewidth',2)


xlabel('time(ms)')
ylabel('Amp')

legend('FD scheme without L-W','traditional FD scheme for L-W','proposed FD scheme 1 for L-W','Time dispersion elimination','proposed FD scheme 2 for L-W')
grid on
axis([ 0 5150 -0.65 1.01])
box on
grid on

% load('figure5bb_2.mat')
%
% % load('figure5c3M8.mat')
%
% % plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
% % xlabel('x/dx')
% % ylabel('Time(ms)')
% % title('')
%
% seis_recordd=seis_record;
% hold on;plot((3:2:nt-2)*1.5/2,seis_recordd((3:2:nt-2),75)./max(abs(seis_recordd((3:2:nt-2),75))),'b','linewidth',2)


% hold on;plot(seis_recorde(:,75)./max(abs(seis_recorde(:,75))),'g','linewidth',2)

% figure;
% % plot(seis_recorda(:,75)./max(abs(seis_recorda(:,75))),'b','linewidth',2);
%
% % hold on;
% plot(1000*(seis_recordb(:,75)./max(abs(seis_recordb(:,75)))-seis_recordc(:,75)./max(abs(seis_recordc(:,75)))),'r','linewidth',2)
% % hold on;plot(seis_recordc(:,75)./max(abs(seis_recordc(:,75)))-(temp./max(abs(temp))),'k','linewidth',2)



% hold on;plot(temp-1*10^-5,'b','linewidth',1);