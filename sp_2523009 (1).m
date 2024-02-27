clc;
clear all;
close all;

Vinf=1;

r=1;
alfa = 0;    
n=4; % number of panels (You can change that number)
step = 2*pi/n;
k = 0 ;

for theta=pi+pi/n:-step:-pi+pi/n
    k = k + 1;
    if r >=0
        X(k)=r*cos(theta);
        Y(k)=r*sin(theta);
    else
        r = -r;
        X(k)=r*cos(theta);
        Y(k)=r*sin(theta);
    end
end


for i=1:1:n

    phi(i)=-alfa+atan2((Y(i+1)-Y(i)),(X(i+1)-X(i)));



    beta(i)=phi(i)+pi/2;

    midx(i)=(X(i+1)+X(i))/2;

    midy(i)=(Y(i+1)+Y(i))/2;

    S(i)=sqrt((Y(i+1)-Y(i))^2+(X(i+1)-X(i))^2);

end


for panel=1:1:n

    others(:,panel)=[1:panel-1 panel+1:n];

    xi=midx(panel);

    yi=midy(panel);

    for i=1:1:n-1

        m=others(i,panel);

        A=-(xi-X(m))*cos(phi(m))-(yi-Y(m))*sin(phi(m));

        B=(xi-X(m))^2+(yi-Y(m))^2;

        C = sin(phi(panel)) * cos(phi(m)) - cos(phi(panel)) * sin(phi(m));

        D=(yi-Y(m))*cos(phi(panel))-(xi-X(m))*sin(phi(panel));

        I(panel,m)=C/2*log((S(m)^2+2*A*S(m)+B)/B)+(D-A*C)/sqrt(B-A^2)*(atan2((S(m)+A),sqrt(B-A^2))-atan2(A,sqrt(B-A^2)));

        J(panel,m)=(D-A*C)/2/sqrt(B-A^2)*log((S(m)^2+2*A*S(m)+B)/B)-C*(atan2((S(m)+A),sqrt(B-A^2))-atan2(A,sqrt(B-A^2)));

    end

    w(panel,1)=Vinf*cos(beta(panel));

end

M=I/2/pi+eye(n)/2;

lambda=-inv(M)*w;

V=Vinf*sin(beta)+lambda'/2/pi*J';

Cp=1-(V/Vinf).^2;

angles=min(beta):0.01:max(beta);

Cp_exact=1-4*sin(angles).^2;




subplot(1,2,1);
    
plot(r*cos(0:0.01:2*pi),r*sin(0:0.01:2*pi),'b',X,Y,'r',midx,midy,'*');
    
legend('Circular Cylinder','Source Panel Method','Control Points')
    
subplot(1,2,2);
    
plot(angles,Cp_exact,'b',beta,Cp,'*');
legend('C_p (analytic)', 'C_p (Source Panel Method)');
title('Cp vs Angle')
ylabel('Cp')
xlabel('Angle')


