function [W1,W2,Kw,Ku,Kd,ww]=gategains(Qu,G1B,H,yu,Qd,G2B,yd)
Cd=0.6;
% Upstream Gate parameters
W1=Qu/(Cd*G1B*sqrt(2*9.81*(H-yu)));
Kw=Cd*G1B*sqrt(2*9.81*(H-yu));
Ku=(-Cd*G1B*W1*sqrt(2*9.81))/(2*sqrt(H-yu));
% Downstream Gate parameters
W2=Qd/(Cd*G2B*sqrt(2*9.81*yd));
Kd=(Cd*G2B*W2*sqrt(2*9.81))/(2*sqrt(yd));

ww=[Qu,W1,Ku,Kw,Qd,W2,Kd];
end
