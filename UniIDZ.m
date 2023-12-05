function [P11,P12,P21,P22,Qm,xx]=UniIDZ(m,n,Sb,L,B,y)

g=9.81;
A=(B+m*y)*y;
P=B+2*(y*sqrt(1+m*m));
dP=(2*sqrt(1+m*m));
R=A/P;
T=B+2*(m*y);
dT=0;
Qm=((R^(2/3))*A*sqrt(Sb)/n)
k=(7/3)-((4*A)/(3*T*P))*dP;
Sf=((Qm^2)*(n^2))/((A^2)*(R^(4/3)));
V=Qm/A;
C=sqrt((9.81*A)/T);
alpha=C+V;
beta=C-V;
gamma=g*(1+k)*Sb;
delta=2*g*Sb/V;

td= L/(alpha);
tu= L/(beta);
                                                                                                                                                                                                                                                                                                                                                                                                                                                              
o=(gamma*L)/(alpha*beta);
r1 = ((alpha*delta)-gamma)/(alpha*(alpha+beta));
r2 = ((beta*delta)-gamma)/(beta*(alpha+beta));
Au = ((alpha*beta*T)/gamma)*(exp(o)-1);
Ad = ((alpha*beta*T)/gamma)*(1-exp(-o));

DenoT=1+(exp(-2*(r1+r2)*L));
NumTbu=(1+(alpha^2/beta^2)*exp(-2*(r1+r2)*L));
NumTbd=(1+(beta^2/alpha^2)*exp(-2*(r1+r2)*L));

bu=(1/(T*alpha))*sqrt((NumTbu/DenoT));
but=((alpha+beta)/(T*alpha*beta))*((exp(-r2*L))/sqrt(DenoT))
bdt=((alpha+beta)/(T*alpha*beta))*((exp(-r1*L))/sqrt(DenoT))
bd=(1/(T*beta))*sqrt((NumTbd/DenoT))

s=tf('s');

% Transfer functions
P11=(1/(Au*s))+bu;
P12=-((1/(Au*s))+but)*exp(-tu*s);
P21=((1/(Ad*s))+bdt)*exp(-td*s);
P22=-((1/(Ad*s))+bd);

xx=[tu,td,Au,Ad,bu,but,bdt,bd];
end
