% 历时 108.221761 秒。
% 历时 114.678604 秒。
% 历时 93.262559 秒。
% 历时 90.642419 秒。
clear
clc %%%%%%%
close all
nt=ceil(3539*(15/28));
eps=.6;     % stability
isnap=60;    % snapshot sampling

load('vv.mat')
c1=flipud(c);
c=c1;

[nz,nx]=size(c1);
v=c1;


nz=nz+50;

vv=zeros(nz,nx);
for ii=1:nz-50
    for jj=1:nx
        vv(ii+50,jj)=v(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vv(ii,jj)=vv(51,jj);
    end
end



clear v
v=vv;



p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=15;  %calculate space increment
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0028; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=51;
xs=nx/2;

h=dx;
r=v*dt/h;
for loop=1:1
    tic

    p=zeros([nz nx]); Vx=p; Vz=p;


    % coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
    % coeff=
    coeff=[1.55417, -0.271537, 0.0738815, -0.0208512, 0.00510231, -0.000914956, 0.000087595];
    % coeff=[1.57242, -0.286902, 0.0846375, -0.0269306, 0.00772107, -0.00167897, 0.000199647]
    % coeff=[1.58345, -0.296503, 0.0918184, -0.031409, 0.00991343, -0.00242128, 0.000326626]


    coeff2=[2.42627, -0.749557, 0.236514, -0.068554, 0.016321, -0.00275895, 0.000243906];
    % coeff2=[2.53363, -0.839206, 0.298263, -0.102656, 0.0305981, -0.00679565, 0.000817839]
    % coeff2=[2.58111, -0.880516, 0.329145, -0.121902, 0.040011, -0.00997917, 0.00136179]


    coeff1=  [ 1.83185, -0.350241, 0.097898, -0.0276368, 0.00670127, -0.00119173, 0.000114199];
    temp1= -3.11498;

    % coeff1= [ 1.86595, -0.37716, 0.11576, -0.0373161, 0.0107582, -0.00236942, 0.000291756]
    % 
    % temp1=-3.15183;

    for it=1:nt-2,

        d1px11=Vx-circshift(Vx,[0 1]);
        d1px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
        d1px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
        d1px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
        d1px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
        d1px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
        d1px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));


        d1pz11=Vz-circshift(Vz,[1 0]);
        d1pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
        d1pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
        d1pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
        d1pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
        d1pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
        d1pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));



        d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
            +coeff(7)*d1px17;
        d1pz=coeff(1)*d1pz11+coeff(2)*d1pz12+coeff(3)*d1pz13+coeff(4)*d1pz14+coeff(5)*d1pz15+coeff(6)*d1pz16...
            +coeff(7)*d1pz17;

        d1px_L=coeff2(1)*d1px11+coeff2(2)*d1px12+coeff2(3)*d1px13+coeff2(4)*d1px14+coeff2(5)*d1px15+coeff2(6)*d1px16...
            +coeff2(7)*d1px17;
        d1pz_L=coeff2(1)*d1pz11+coeff2(2)*d1pz12+coeff2(3)*d1pz13+coeff2(4)*d1pz14+coeff2(5)*d1pz15+coeff2(6)*d1pz16...
            +coeff2(7)*d1pz17;

        Q1=-v.^2.*(d1px_L+d1pz_L);


        d2px=(circshift(Q1,[0 -1])+circshift(Q1,[0 1]));
        d2pz=(circshift(Q1,[-1 0])+circshift(Q1,[1 0]));

        d2pxz=v.^2.*(d2px+d2pz-4*Q1);


        p=p-dt*v.^2.*(d1px+d1pz)./h +d2pxz*dt^3/h^3/24;

        p(zs,xs)= p(zs,xs)+src(it)*dt^2;

        [p,p]=spongeABC(p,p,nx,nz,50,50,0.007);

        d1px=(circshift(p,[0 -1])-p);


        d1pz=(circshift(p,[-1])-p);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d2pz11=(circshift(p,[-1 0])+circshift(p,[1 0]));
        d2pz12=(circshift(p,[-2 0])+circshift(p,[2 0]));
        d2pz13=(circshift(p,[-3 0])+circshift(p,[3 0]));
        d2pz14=(circshift(p,[-4 0])+circshift(p,[4 0]));
        d2pz15=(circshift(p,[-5 0])+circshift(p,[5 0]));
        d2pz16=(circshift(p,[-6 0])+circshift(p,[6 0]));
        d2pz17=(circshift(p,[-7 0])+circshift(p,[7 0]));


        d2pz=temp1*p+coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
            +coeff1(7)*d2pz17;




        d2px11=(circshift(p,[0 -1])+circshift(p,[0 1]));
        d2px12=(circshift(p,[0 -2])+circshift(p,[0 2]));
        d2px13=(circshift(p,[0 -3])+circshift(p,[0 3]));
        d2px14=(circshift(p,[0 -4])+circshift(p,[0 4]));
        d2px15=(circshift(p,[0 -5])+circshift(p,[0 5]));
        d2px16=(circshift(p,[0 -6])+circshift(p,[0 6]));
        d2px17=(circshift(p,[0 -7])+circshift(p,[0 7]));
        %
        %
        d2px=temp1*p+coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
            +coeff1(7)*d2px17;


        d2pxz=v.^2.*(d2px+d2pz);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d3pz=circshift(d2pxz,[-1])-d2pxz;
        Vz=Vz-dt*d1pz/h -(d3pz)*dt^3/24/h^3;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d3px=circshift(d2pxz,[0 -1])-d2pxz;
        Vx=Vx-dt*d1px/h -(d3px)*dt^3/24/h^3;
        %%%%%%%%%%%%%%%%%%%%%%%
        [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.007);

        if rem(it,isnap)== 0,
            imagesc(x,z,p), axis equal
            colormap gray
            xlabel('x'),ylabel('z')
            title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
            drawnow
        end

        seis_record(it,:)=p(zs,:);
    end

end
toc
save('figure8b.mat')

