

clear;clc;close all
M=7

f=  [ 1.83185, -0.350241, 0.097898, -0.0276368, 0.00670127, -0.00119173, 0.000114199];
f1= -3.11498;

  c=[1.55417, -0.271537, 0.0738815, -0.0208512, 0.00510231, -0.000914956, 0.000087595];


    d=[2.42627, -0.749557, 0.236514, -0.068554, 0.016321, -0.00275895, 0.000243906];



h=10;
v=2500;
v=5000;
% v=8500;
tau=0.001
dt=tau;
r=v*tau/h;
kh=linspace((pi)/(100),(pi),100);

for jj=1:5
    xita=(jj-1)*pi/16;
    d1kz=0;
    for m=1:M
        d1kz=d1kz+2*c(m)*sin((m-0.5)*kh*sin(xita));
    end
    %     d1kz=d1kz./h;

    d1kx=0;
    for m=1:M
        d1kx=d1kx+2*c(m)*sin((m-0.5)*kh*cos(xita));
    end
    %     d1kx=d1kx./h;

    d1kz_1=0;
    for m=1:M
        d1kz_1=d1kz_1+2*d(m)*sin((m-0.5)*kh*sin(xita));
    end
    %     d1kz=d1kz./h;

    d1kx_1=0;
    for m=1:M
        d1kx_1=d1kx_1+2*d(m)*sin((m-0.5)*kh*cos(xita));
    end
    %     d1kx=d1kx./h;


    d2kx=f1;
    for m=1:M
        d2kx=d2kx+2*f(m)*(cos(m*kh*cos(xita))    );
    end
    %     d2kx=d2kx./h^2;

    d2kz=f1;
    for m=1:M
        d2kz=d2kz+2*f(m)*(cos(m*kh*sin(xita))   );
    end


    %      d2kz=d2kz./h^2;
    %%%%%%%%%%%%%%%%%%%
    temp1=d2kx+d2kz;
    temp2=-4+2*cos(kh*cos(xita))+2*cos(kh*sin(xita));

    g11=-d1kx.*(2*sin(0.5*kh*cos(xita))+r^2*1/24*2*sin(0.5*kh*cos(xita)).*temp1);
    g12=-d1kz.*(2*sin(0.5*kh*sin(xita))+r^2*1/24*2*sin(0.5*kh*sin(xita)).*temp1);
    %     g1=v^2*dt^2*(g11+g12);

    g1=r^2*(g11+g12);

    g11=-d1kx_1.*(2*sin(0.5*kh*cos(xita))+r^2*1/24*2*sin(0.5*kh*cos(xita)).*temp1);
    g12=-d1kz_1.*(2*sin(0.5*kh*sin(xita))+r^2*1/24*2*sin(0.5*kh*sin(xita)).*temp1);


    g1_1=r^2*(g11+g12);


    atemp=(g1)+r^2/24*(g1_1.*temp2);


    %%%%%%%%%%%%%%%


    a1=real(h/v*(1./(acos(1+1/2*atemp)./(kh*r))-1));

    if (jj==1)
        figure;plot(100*kh/(pi),a1,'m','linewidth',2)
        hold on
    elseif jj==2
        plot(100*kh/(pi),a1,'r--','linewidth',2)
    elseif jj==3
        plot(100*kh/(pi),a1,'c:','linewidth',2)
    elseif jj==4
        plot(100*kh/(pi),a1,'k-.','linewidth',2)
    else
        plot(100*kh/(pi),a1,'b','linewidth',2)
    end

end

axis([0 100 -0.5*10^-5 7*10^-5])
% axis([0 100 -0.5*10^-5 3*10^-5])
grid on
xlabel('percentage of kh')
legend('\theta=0', '\theta=π/16','\theta=2π/16','\theta=3π/16','\theta=4π/16')

ylabel('\xi_1_ (\theta)')