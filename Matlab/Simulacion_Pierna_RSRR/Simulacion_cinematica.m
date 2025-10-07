close all
clc
clear
num=100;
Tt=[];
cir=linspace(0,pi,num);
cir_int=linspace(pi/12,pi-pi/12,num);
Rim=linspace(0,pi/6,10);
Rim_int=linspace(pi+pi/4,2*pi-pi/4,50);
alpha_bpt=[pi/12 2*pi/9 13*pi/36 pi/2];
alpha_pt=[pi/12 pi/6 3*pi/12];
lp2ct=[42.6 50 57.4];
posr=zeros(2,2,7);
rt=75;
rt_int=60;
rc_ext=40;
rc_int=20;
rc_int2=35;
raxle=15;
r1=52.5;
r2=27.5;
lp2=65;

%Configuracion de la rueda
i=1;
j=1;
k=2;

posr(:,:,1)=[cos(pi/12) -sin(pi/12);
            sin(pi/12) cos(pi/12)];
posr(:,:,2)=[cos(75*pi/180) -sin(75*pi/180);
            sin(75*pi/180) cos(75*pi/180)];
posr(:,:,3)=[cos(135*pi/180) -sin(135*pi/180);
            sin(135*pi/180) cos(135*pi/180)];
posr(:,:,4)=[cos(35*pi/180) -sin(35*pi/180);
            sin(35*pi/180) cos(35*pi/180)];
posr(:,:,5)=[cos(pi+pi/4) -sin(pi+pi/4);
            sin(pi+pi/4) cos(pi+pi/4)];
posr(:,:,6)=[cos(pi+105*pi/180) -sin(pi+105*pi/180);
            sin(pi+105*pi/180) cos(pi+105*pi/180)];
posr(:,:,7)=[cos(pi+165*pi/180) -sin(pi+165*pi/180);
            sin(pi+165*pi/180) cos(pi+165*pi/180)];

alpha_bp=alpha_bpt(i);
alpha_p=alpha_pt(j);

mqp1=[31.457 43.714 39.862];%Proximal_11
mqp2=[61.457 58.714 54.862];%Proximal_21
%Archivo que contiene la trayectoria
load('Walk_50.mat')
%xy=[-50 0 50;-100 -75 -100];
%mov = VideoWriter('Cinematic');
%while true
N=length(xy);%100;
alpha21=[0;0];
alpha22=[0;0];
alpha23=[0;0];
%N=length(xy);
%for k=1:3
 lp2c=lp2ct(k);
 Ml=sqrt(lp2^2-(lp2c-23.75)^2);
 dqp=mqp1(k)/N;
 dqr=2*pi/N;
 qr=0;
 qp=0;%mqp1(2);
 
 T=2;
NP=100;
fnum=1;
for i=1:NP
%% for idx=1:N

     
theta_p1=pi+ pi/8;
theta_p2=2*pi-pi/8;
%  
%  
eCx=0;
eCy=-100;
eA=40;
eB=20;

%  th_bar=pi-pi*(i-1)/(NP-1);
%  x_t(i)=eCx + eA*cos(th_bar);
%  y_t(i)=eCy + eB*sin(th_bar);
%   
tiempo(i)=T*(i-1)/(NP-1);
if i<(NP/2)
    
    CLOID=(i/(NP/2)-(1/(2*pi))*sin(2*pi*i/(NP/2)));
    th_bar=pi-pi*CLOID;
    x_t(i)=eCx + eA*cos(th_bar);
    y_t(i)=eCy + eB*sin(th_bar);

else
    CLOID=((i-NP/2)/(NP/2)-(1/(2*pi))*sin(2*pi*(i-NP/2)/(NP/2)));
    x_t(i)=(eCx+eA) - 2*eA*CLOID ;
    y_t(i)=eCy;

end
 
%for idx=1:5:N%length(xy(1,:))
%Obtencion de las variables articulares con las Redes
qp=Fcn_Red_Qp_1([x_t(i), y_t(i)]');   
qr=Fcn_Red_Qr_1([x_t(i), y_t(i)]');

qp_h(i)=qp;
qr_h(i)=qr;


%qp=qp+dqp;
%qr=qr+dqr;
[lp1, lp12, h1, h2, lp2t, alpha_2, alpha_3, alpha_4, rot, r11, c1, b1, p2c, Ti, Ti2, ut, a_p]=Cin_dir_2(qr, qp, rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tt=[Tt Ti];

RuedaH_1=rot*[rt*cos(cir) rt_int*cos(flip(cir));rt*sin(cir) rt_int*sin(flip(cir))];
RuedaH_1=[RuedaH_1;zeros(1,length(RuedaH_1))];
Rim_1=rot*posr(:,:,1)*[rt_int*cos(Rim),rc_int*cos(flip(Rim));rt_int*sin(Rim),rc_int*sin(flip(Rim))];
Rim_1=[Rim_1;zeros(1,length(Rim_1))];
Rim_2=rot*posr(:,:,2)*[rt_int*cos(Rim),rc_int*cos(flip(Rim));rt_int*sin(Rim),rc_int*sin(flip(Rim))];
Rim_2=[Rim_2;zeros(1,length(Rim_2))];
Rim_3=rot*posr(:,:,3)*[rt_int*cos(Rim),rc_int*cos(flip(Rim));rt_int*sin(Rim),rc_int*sin(flip(Rim))];
Rim_3=[Rim_3;zeros(1,length(Rim_3))];
Rim_4=rot*posr(:,:,4)*[rc_int*cos(cir) rc_int2*cos(Rim_int);rc_int*sin(cir) rc_int2*sin(Rim_int)];
Rim_4=[Rim_4;zeros(1,length(Rim_4))];
Eje=[raxle*cos(cir) raxle*cos(pi+cir);raxle*sin(cir) raxle*sin(pi+cir)];
Eje=[Eje;zeros(1,length(Eje))];


RuedaH_2=p2c+rot*[rt*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir) rt_int*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(cir));-rt*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir) -rt_int*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(cir))];
RuedaH_2=[RuedaH_2;zeros(1,length(RuedaH_2))];
Rim2_1=p2c+rot*posr(:,:,5)*[rt_int*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),rc_ext*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(Rim));-rt_int*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),-rc_ext*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(Rim))];
Rim2_1=[Rim2_1;zeros(1,length(Rim2_1))];
Rim2_2=p2c+rot*posr(:,:,6)*[rt_int*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),rc_ext*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(Rim));-rt_int*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),-rc_ext*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+flip(Rim))];
Rim2_2=[Rim2_2;zeros(1,length(Rim2_2))];
Rim2_3=p2c+rot*posr(:,:,7)*[rt_int*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),rc_ext*cos(flip(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim));-rt_int*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim),-rc_ext*sin(flip(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+Rim))];
Rim2_3=[Rim2_3;zeros(1,length(Rim2_3))];
Rim2_4=p2c+rot*[rc_ext*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir_int),rc_int2*cos(flip(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir_int));-rc_ext*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir_int),-rc_int2*sin(flip(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir_int))];
Rim2_4=[Rim2_4;zeros(1,length(Rim2_4))];



SemC_1=rot*[rt*cos(cir);rt*sin(cir)];
SemC_2=p2c+rot*[rt*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir);-rt*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir)];
Cbp=rot*[r2*cos(2*cir);r2*sin(2*cir)];
Cp1=p2c+rot*[r1*cos(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir);-r1*sin(pi/2+alpha_4+alpha_3+alpha_2+alpha_bp+cir)];

fig=figure(fnum);
%clf
hold on
%Parte fija de la rueda
ColorP= [0.8157,    0.8235,    0.1176];
ColorP= [0.4, 0.4, 0.4];
ColorP=ColorP+(min(1-ColorP)-0.00)*(NP-i)/NP; 

fill3(RuedaH_1(1,:),RuedaH_1(2,:),RuedaH_1(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim_1(1,:),Rim_1(2,:),Rim_1(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim_2(1,:),Rim_2(2,:),Rim_2(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim_3(1,:),Rim_3(2,:),Rim_3(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim_4(1,:),Rim_4(2,:),Rim_4(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
%fill3(Rim_4(1,:),Rim_4(2,:),Rim_4(3,:),[180/255,180/255,180/255],'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Eje(1,:),Eje(2,:),Eje(3,:),'w','EdgeColor','none')%, FaceAlpha=0.2)



% ColorP= [0.8157,    0.8235,    0.1176];
% 
% ColorP=ColorP+(min(1-ColorP)-0.00)*(NP-i)/NP; 


%Parte movil de la rueda
fill3(RuedaH_2(1,:),RuedaH_2(2,:),RuedaH_2(3,:), ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim2_1(1,:),Rim2_1(2,:),Rim2_1(3,:), ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim2_2(1,:),Rim2_2(2,:),Rim2_2(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim2_3(1,:),Rim2_3(2,:),Rim2_3(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
fill3(Rim2_4(1,:),Rim2_4(2,:),Rim2_4(3,:),ColorP,'EdgeColor','none')%, FaceAlpha=0.5)
%Proximal 1

neg=0.4
plot3([r11(1) c1(1)],[r11(2) c1(2)], [0 0], 'LineWidth', 4, 'Color', ColorP-neg)
plot3([r11(1) c1(1)],[r11(2) c1(2)], [0 0], 'o', 'LineWidth', 1, 'Color',  ColorP-neg)
%Proximal 2
plot3([0 0], [0 0], [0 Ml-qp], 'LineWidth', 4, 'Color',  ColorP-neg)
plot3([0 23.75*sin(pi-qr)], [0 23.75*cos(pi-qr)], [Ml-qp Ml-qp], 'LineWidth', 2, 'Color',  ColorP-neg)
plot3([23.75*sin(pi-qr) b1(1)],[23.75*cos(pi-qr) b1(2)],[Ml-qp 0], 'LineWidth', 2, 'Color',  ColorP-neg)
plot3(b1(1),b1(2),0,'o', 'LineWidth', 1, 'Color',  ColorP-neg)
%Trayectoria
plot3(Tt(1,:),Tt(2,:),0*Tt(1,:))
%title(sprintf('x=%f  qr=%f\ny=%f qp=%f',Ti(1),qr,Ti(2),qp))

grid on
%view(3)%(180,-90)
pbaspect([1 1 1])
xlim([-150 150])
ylim([-150 150])
zlim([-1 65])
%pause(0.5)


%  F = getframe(fig);
%  open(mov);
%  writeVideo(mov,F);
end
fnum=fnum+1;

save trajectory Tt


set(fig,'Color',[1 1 1])
TamLet=14;
xlabel('$x$ [mm]','Interpreter','latex','FontSize',TamLet)
ylabel('$y$ [mm]','Interpreter','latex','FontSize',TamLet)



fig2=figure(fnum);
plot(tiempo, x_t, tiempo, y_t,'--')
axis square 
grid on
xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
ylabel('[mm]','Interpreter','latex','FontSize',TamLet)

legend({'$t_x$','$t_y$'},'Interpreter','latex', 'FontSize',TamLet)

set(fig2,'Color',[1 1 1])
fnum=fnum+1;


% fig3=figure(fnum);
% plot(tiempo, qp_h)
% %axis([0 2 -45 10])
% axis square 
% grid on
% xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
% ylabel('$q_p$ [mm]','Interpreter','latex','FontSize',TamLet)
% set(fig3,'Color',[1 1 1])
% fnum=fnum+1;
% 
% fig4=figure(fnum);
% plot(tiempo, qr_h*180/pi)
% %axis([0 2 -45 10])
% axis square 
% grid on
% xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
% ylabel('$q_r$ [deg]','Interpreter','latex','FontSize',TamLet)
% set(fig4,'Color',[1 1 1])
% fnum=fnum+1;


fig5=figure(fnum);
plot(tiempo(2:NP), diff(x_t)/(T/NP), tiempo(2:NP), diff(y_t)/(T/NP),'--')
%axis([0 2 -45 10])
axis square 
grid on
xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
ylabel('[mm/s]','Interpreter','latex','FontSize',TamLet)

legend({'$\dot{t}_x$','$\dot{t}_y$'},'Interpreter','latex', 'FontSize',TamLet)
set(fig5,'Color',[1 1 1])
fnum=fnum+1;
% 
% 
% fig6=figure(fnum);
% plot(tiempo(2:NP), diff(qp_h)/(T/NP))
% axis square 
% %axis([0 2 -45 10])
% grid on
% xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
% ylabel('$\dot{q}_p$ [mm/s]','Interpreter','latex','FontSize',TamLet)
% set(fig6,'Color',[1 1 1])
% fnum=fnum+1;

% fig7=figure(fnum);
% plot(tiempo(2:NP), diff(qr_h)/(T/NP))
% %axis([0 2 -45 10])
% axis square 
% grid on
% xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
% ylabel('$\dot{q}_r$ [rad/s]','Interpreter','latex','FontSize',TamLet)
% set(fig7,'Color',[1 1 1])
% fnum=fnum+1;

fig8=figure(fnum);
plot(tiempo, qr_h, tiempo(2:NP), diff(qr_h)/(T/NP),'-.g')
axis square 
grid on
xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
%ylabel('[mm]','Interpreter','latex','FontSize',TamLet)

legend({'$q_r$ [rad]','$\dot{q}_r$ [rad/s]'}, 'Interpreter','latex','FontSize',12)

set(fig8,'Color',[1 1 1])
fnum=fnum+1;

fig9=figure(fnum);
plot(tiempo, qp_h, tiempo(2:NP), diff(qp_h)/(T/NP),'-.g')
axis square 
grid on
xlabel('t [s]','Interpreter','latex','FontSize',TamLet)
%ylabel('[mm]','Interpreter','latex','FontSize',TamLet)

legend({'$q_p$ [mm]','$\dot{q}_p$ [mm/s]'},'Interpreter','latex', 'FontSize',TamLet)

set(fig9,'Color',[1 1 1])





