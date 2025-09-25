function [lp1 lp12 h1 h2 lp2t alpha_2 alpha_3 alpha_4 rot r11 c1 b1 p2c Ti Ti2 ut a_p]=Cin_dir_2(qr, qp, rt, r1, r2, lp2c, lp2, alpha_bp, alpha_p, Ml)
lp1=sqrt(r1^2+r2^2-2*r1*r2*cos(alpha_bp-alpha_p));
lp12=sqrt(r1^2+lp2c^2-2*r1*lp2c*cos(pi/2-alpha_p));

h1=sqrt(lp2^2-(Ml-qp)^2)+23.75;
h2=sqrt(r2^2+h1^2-2*r2*h1*cos(pi/2-alpha_bp));
lp1t=sqrt(r1^2+rt^2-2*r1*rt*cos(pi-alpha_p));
lp2t=sqrt(rt^2+lp2c^2);
alpha_2=pi-(acos(-(h1^2-h2^2-r2^2)/(2*r2*h2))+acos(-(lp12^2-lp1^2-h2^2)/(2*lp1*h2)));
alpha_3=pi-acos(-(h2^2-lp1^2-lp12^2)/(2*lp1*lp12));
alpha_4=pi-acos(-(r1^2-lp12^2-lp2c^2)/(2*lp12*lp2c));
ut=acos(-(lp2t^2-lp1t^2-lp12^2)/(2*lp1t*lp12));

alpha_5=alpha_3+ut;

rot=[cos(qr) -sin(qr);
    sin(qr) cos(qr)];


r11=rot*[r2*cos(alpha_bp);-r2*sin(alpha_bp)];
c1=r11+rot*[lp1*cos(alpha_2+alpha_bp);-lp1*sin(alpha_2+alpha_bp)];
b1=c1+rot*[lp12*cos(alpha_3+alpha_2+alpha_bp);-lp12*sin(alpha_3+alpha_2+alpha_bp)];
p2c=b1+rot*[lp2c*cos(alpha_4+alpha_3+alpha_2+alpha_bp);-lp2c*sin(alpha_4+alpha_3+alpha_2+alpha_bp)];
Ti=p2c+rot*[rt*cos(-pi/2+alpha_4+alpha_3+alpha_2+alpha_bp);-rt*sin(-pi/2+alpha_4+alpha_3+alpha_2+alpha_bp)];
Ti2=c1+rot*[lp1t*cos(alpha_5+alpha_2+alpha_bp);-lp1t*sin(alpha_5+alpha_2+alpha_bp)];
a_p=[lp1*cos(alpha_2+alpha_bp);-lp1*sin(alpha_2+alpha_bp)];
end
