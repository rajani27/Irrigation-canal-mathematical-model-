function [P11,P12,P21,P22,xx]=IDZ(m,n,Sb,L,y2,B,y)

A=(B+m*y)*y;
P=B+2*(y*sqrt(1+m*m));
dP=(2*sqrt(1+m*m))
R=A/P;
T=B+2*(m*y)
dT=0;
Qm=((R^(2/3))*A*sqrt(Sb)/n)
k=(7/3)-((4*A)/(3*T*P))*dP;
Sf=((Qm^2)*(n^2))/((A^2)*(R^(4/3)));
V=Qm/A;
C=sqrt((9.81*A)/T);
al=C+V;
be=C-V;
F=sqrt((V^2)*T/(9.81*A));
Sx=(Sb-Sf)/(1-(F^2));
d=(2*9.81/V)*(Sf-(Sx*(F^2)));
g=((C*C/T)*dT)+(9.81*(((1+k)*Sb)-(1+k-(k-2)*F^2)*Sx));

% Backwater effect

A2=(B+m*y2)*y2;
P2=B+2*(y2*sqrt(1+m*m));
R2=A2/P2;
Sf2=((Qm.^2)*(n.^2))/((A2.^2)*(R2.^(4/3)));
T2=B+2*(m*y2)
V2=Qm/A2;
C2=sqrt((9.81*A2)/T2);
F2=sqrt((V2^2)*T2/(9.81*A2));
Sx2=(Sb-Sf2)/(1-(F2^2));
if Sx2~=0;
x1=max(L-((y2-y)/Sx2),0);
else
x1=L;
end
x2=L-x1;
%Equivalent backwater areas
Ad=(1-exp(-g*x1/(T*(C*C-V*V))))*(T^2*((C^2)-(V^2)))/g;  
Au=(exp(g*x1/(T*(C*C-V*V)))-1)*(T^2*((C^2)-(V^2)))/g;  
%and the delays
td=x1/(al);
tu=x1/(be);
t_d=x2/(al);
t_u=x2/(be);
%Equivalent areas are given by
A_d=(1-exp(-g*x2/(T*(C*C-V*V))))*(T^2*((C^2)-(V^2)))/g;
A_u=(exp(g*x2/(T*(C*C-V*V)))-1)*(T^2*((C^2)-(V^2)))/g;
AD=A_d*(1+(Ad/A_u));
AU=Au*(1+(A_u/Ad));
tD=td+t_d;
tU=tu+t_u;

p11=tf(1,[AU 0]);
p22=tf(-1,[AD 0]);
eu=tf([1],[1],'InputDelay',tU);
p12=-p11*eu;
ed=tf(1,1,'InputDelay',tD);
p21=-p22*ed;

%High Frequencies approximation

%for u/s part
af=(T*(2+(k-1)*F*F)*Sb)/(A*F*(1-(F^2)))
g11=(1/(T*C*(1-F)))*sqrt((1+exp(af*x1)*((1-F)/(1+F))^2)/(1+exp(af*x1)));
g12=(2/(T*C*(1-F^2)))*(exp(-g*x1/(2*T*(C*C-V*V)))/sqrt(1+exp(af*x1)))
g21=(2/(T*C*(1-F^2)))*(exp(g*x1/(2*T*(C*C-V*V)))/sqrt(1+exp(af*x1)))
g22=(1/(T*C*(1+F)))*sqrt((1+exp(af*x1)*((1+F)/(1-F))^2)/(1+exp(af*x1)))

%for d/s part
alf=(T/(A*F*(1-(F^2))))*[(2+(k-1)*F^2)*Sb-(2+(k-1)*(F^2)-(A*dT/(T^2)+k-2)*F^4)*Sx]
h11=(1/(T*C*(1-F)))*sqrt((1+exp(alf*x2)*((1-F)/(1+F))^2)/(1+exp(alf*x2)))
h12=(2/(T*C*(1-F^2)))*(exp(-g*x2/(2*T*(C*C-V*V)))/sqrt(1+exp(af*x2)))
h21=(2/(T*C*(1-F^2)))*(exp(g*x2/(2*T*(C*C-V*V)))/sqrt(1+exp(af*x2)))
h22=(1/(T*C*(1+F)))*sqrt((1+exp(alf*x2)*((1+F)/(1-F))^2)/(1+exp(alf*x2)))

T11=g11+((g12*g21)/(h11+g22))
T12=-((g12*h12)/(h11+g22))*eu
T21=((g21*h21)/(h11+g22))*ed
T22=-(h11+((h12*h21)/(h11+g22)))

%global approximate model
P11=p11+T11
P12=p12+T12
P21=(p21+T21)
P22=(p22+T22)

xx=[Qm,tU,tD,AU,AD,x1,x2,T11,T22];
end
