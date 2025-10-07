clear
clc
load('Test.mat')
%Parameters definition
alpha=0.1;
R_c1=27.5;
L_1=65;
L_2=62.4295604118095;
L_3=25;
R_c2=23.75;
aux1=sqrt(52.5^2+75^2-2*75*52.5*cos(11*pi/12));
l_t=sqrt(75^2+50^2);
l_p=59.464;
delta=-5/12*pi;
gamma=-7/12*pi;
DataQ=[];
DataT=[];
imax=200;
%Inverse kinematics and Newton-Raphson
for k=1:length(Qr)
    
    Ti_d=T(:,k);
    qr=Qr(k)-0.1;
    qp=Qp(k)-1;
    
    for idx=1:imax
        rho_11=qr+gamma;
        sigma_1=rho_11+delta;
        theta_1=atan(sqrt((L_1/(l_p-qp))^2-1));
        S_p=R_c2+sqrt(L_1^2-(l_p-qp)^2);
        Q_1=S_p*cos(sigma_1)-R_c1*cos(rho_11);
        Q_2=S_p*sin(sigma_1)-R_c1*sin(rho_11);
        P_1=(S_p^2+R_c1^2-2*S_p*R_c1*cos(sigma_1-rho_11)+L_2^2-L_3^2)/(2*L_2);
        P_2=(S_p^2+R_c1^2-2*S_p*R_c1*cos(sigma_1-rho_11)+L_3^2-L_2^2)/(2*L_3);
        theta_2=2*atan((Q_2-sqrt(Q_1^2+Q_2^2-P_1^2))/(P_1+Q_1));
        theta_3=2*atan((Q_2+sqrt(Q_1^2+Q_2^2-P_2^2))/(P_2+Q_1));
        th=acos(-(-l_t^2+aux1^2-L_2^2)/(2*L_2*l_t))-pi;
        r_11=R_c1*[cos(rho_11+pi/2);sin(rho_11+pi/2)];
        c_1=L_3*[-sin(theta_3);cos(theta_3)];
        b_1=L_2*[-sin(theta_2);cos(theta_2)];
        R_th=[cos(th) -sin(th);
              sin(th) cos(th)];
        K=[-L_2*sin(theta_2) -L_3*sin(theta_3);
            L_2*cos(theta_2) L_3*cos(theta_3)];
        L=[-Q_2 cos(sigma_1)/tan(theta_1);
            Q_1 sin(sigma_1)/tan(theta_1)];
        M=[-l_t*sin(theta_2+th)-L_2*sin(theta_2) -L_3*sin(theta_3);
            l_t*cos(theta_2+th)+L_2*cos(theta_2) L_3*cos(theta_3)];
        N=[-R_c1*sin(rho_11) 0;
            R_c1*cos(rho_11) 0];
        J=M*inv(K)*L+N;
        T_i=((l_t/L_2)*R_th+[1 0;0 1])*b_1+c_1+r_11;
        Dq=inv(J)*(T_i-Ti_d);
     
        qr1=qr-alpha*Dq(1);
        qp1=qp-alpha*Dq(2);
        nE=norm([(qr1-qr),(qp1-qp)]);
        if nE<=0.005
            break
        end
        if ~isreal(qr) || ~isreal(qp)
            break
        end
        qr=qr1;
        qp=qp1;
    end

if isreal(qr) && isreal(qp) && isreal(T_i)
    DataQ=[DataQ [Qr(k);Qp(k);k;idx]];
    DataT=[DataT [T_i Ti_d]];
end
end

Err=[];
for i=1:2:length(DataT)-1
    M1=norm(DataT(:,i));
    M2=norm(DataT(:,i+1));
    Err=[Err norm(M2-M1)];
end
figure
histogram(Err)
title('Error with Newton-Raphson algorithm')
xlabel('Absolute error in mm')
ylabel('Frequency')

