clear
clc
load('Test.mat');

%Parameters definition
rt=75;
r1=52.5;
r2=27.5;
lp2c=50;
lp2=65;
alpha_bp=pi/12;
alpha_p=pi/12;
Ml=sqrt(lp2^2-(lp2c-23.75)^2);
%Inverse kinematics with NN
Qp_r=Fcn_Red_Qp_1(T);
Qr_r=Fcn_Red_Qr_1(T);
%Forward kinematics with Q given by NN
for idx=1:length(Qp_r)
[lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut]=Cin_dir_2(Qr_r(idx), Qp_r(idx), rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml);
Tr(:,idx)=Ti;
end
%Error
T_e=Tr-T;
k=1;
for idx=1:length(T_e)
    Et=sqrt(T_e(1,idx)^2+T_e(2,idx)^2);
    %Filter error greater than 0.1 mm in distance
    if Et>=0.1
        Ea(:,k)=Et;
        k=k+1;
    end
end
