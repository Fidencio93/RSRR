close all
clear
clc
load('Test.mat');


Threshold=0.1;
rt=75;
r1=52.5;
r2=27.5;
lp2c=50;
lp2=65;
alpha_bp=pi/12;
alpha_p=pi/12;
Ml=sqrt(lp2^2-(lp2c-23.75)^2);
Qp_r1=Fcn_Red_Qp_1(T);
Qr_r1=Fcn_Red_Qr_1(T);
for idx=1:length(Qp_r1)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r1(idx), Qp_r1(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr1(:,idx)=Ti;
end
Qp_r2=Fcn_Red_Qp_2(T);
Qr_r2=Fcn_Red_Qr_2(T);
for idx=1:length(Qp_r2)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r2(idx), Qp_r2(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr2(:,idx)=Ti;
end
Qp_r3=Fcn_Red_Qp_3(T);
Qr_r3=Fcn_Red_Qr_3(T);
for idx=1:length(Qp_r3)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r3(idx), Qp_r3(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr3(:,idx)=Ti;
end
Qp_r4=Fcn_Red_Qp_4(T);
Qr_r4=Fcn_Red_Qr_4(T);
for idx=1:length(Qp_r4)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r4(idx), Qp_r4(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr4(:,idx)=Ti;
end
Qp_r5=Fcn_Red_Qp_5(T);
Qr_r5=Fcn_Red_Qr_5(T);
for idx=1:length(Qp_r5)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r5(idx), Qp_r5(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr5(:,idx)=Ti;
end
Qp_r6=Fcn_Red_Qp_6(T);
Qr_r6=Fcn_Red_Qr_6(T);
for idx=1:length(Qp_r6)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r6(idx), Qp_r6(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr6(:,idx)=Ti;
end
T_e1=Tr1-T;
k=1;
for idx=1:length(T_e1)
    Et1=sqrt(T_e1(1,idx)^2+T_e1(2,idx)^2);
    if Et1>=Threshold
        Ea1(:,k)=Et1;%T_e(:,idx);
        k=k+1;
    end
end
T_e2=Tr2-T;
k=1;
for idx=1:length(T_e2)
    Et2=sqrt(T_e2(1,idx)^2+T_e2(2,idx)^2);
    if Et2>=Threshold
        Ea2(:,k)=Et2;%T_e(:,idx);
        k=k+1;
    end
end
T_e3=Tr3-T;
k=1;
for idx=1:length(T_e3)
    Et3=sqrt(T_e3(1,idx)^2+T_e3(2,idx)^2);
    if Et3>=Threshold
        Ea3(:,k)=Et3;%T_e(:,idx);
        k=k+1;
    end
end
T_e4=Tr4-T;
k=1;
for idx=1:length(T_e4)
    Et4=sqrt(T_e4(1,idx)^2+T_e4(2,idx)^2);
    if Et4>=Threshold
        Ea4(:,k)=Et4;%T_e(:,idx);
        k=k+1;
    end
end
T_e5=Tr5-T;
k=1;
for idx=1:length(T_e5)
    Et5=sqrt(T_e5(1,idx)^2+T_e5(2,idx)^2);
    if Et5>=Threshold
        Ea5(:,k)=Et5;%T_e(:,idx);
        k=k+1;
    end
end
T_e6=Tr6-T;
k=1;
for idx=1:length(T_e6)
    Et6=sqrt(T_e6(1,idx)^2+T_e6(2,idx)^2);
    if Et6>=Threshold
        Ea6(:,k)=Et6;%T_e(:,idx);
        k=k+1;
    end
end
cols=[0.1 0.1:0.1:1 1];
figure(1)
histogram(Ea1,cols)
title('Error from Net 1')
xlabel('Absolute error in mm')
ylabel('Frequency')

figure(2)
histogram(Ea2,cols)
title('Error from Net 2')
xlabel('Absolute error in mm')
ylabel('Frequency')

figure(3)
histogram(Ea3,cols)
title('Error from Net 3')
xlabel('Absolute error in mm')
ylabel('Frequency')

figure(4)
histogram(Ea4,cols)
title('Error from Net 4')
xlabel('Absolute error in mm')
ylabel('Frequency')

figure(5)
histogram(Ea5,cols)
title('Error from Net 5')
xlabel('Absolute error in mm')
ylabel('Frequency')

figure(6)
histogram(Ea6,cols)
title('Error from Net 6')
xlabel('Absolute error in mm')
ylabel('Frequency')
