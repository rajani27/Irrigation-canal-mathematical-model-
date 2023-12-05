function [P11,P12,P21,P22,xx]=UniSV(m,n,Sb,L,B,y)

g=9.81;
A=(B+m*y)*y;
P=B+2*(y*sqrt(1+m*m));
dP=(2*sqrt(1+m*m));
R=A/P;
T=B+2*(m*y);
dT=0;
Qm=((R^(2/3))*A*sqrt(Sb)/n);
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


o=gamma*L/(alpha*beta);
r=(alpha-beta)/(alpha*beta);
a=delta/gamma;
b=(T*(exp(gamma*L/(alpha*beta))-1)^2);

s=tf('s');

B11=(1/b)*(r*(exp(o)*(1-o)-1) + a*(exp(o)*(exp(o)-1)+exp(o)*(1-2*o)-1));
B12=(1/b)*(r*(o*exp(o)+(1-exp(o)))+a*(o*(1+exp(o))+2*(1-exp(o))));
B21=(exp(o)/b)*(r*(exp(o)-1-o) + a*((exp(o)-1)+exp(o)*(1-o)-1-o));
B22=(1/b)*(r*exp(o)*(1-exp(o)+o) + a*(1+2*o*exp(o)-exp(o)*exp(o)));

a011=(gamma/(alpha*beta*T))*(1/(exp(o)-1));
a012=-(gamma/(alpha*beta*T))*(1/(exp(o)-1));
a021=(gamma/(alpha*beta*T))*(exp(o)/(exp(o)-1));
a022=-(gamma/(alpha*beta*T))*(exp(o)/(exp(o)-1));


b11= B11;
b12= B12 + tu*a012;
b21= B21 + td*a021;
b22= B22;

t11=b11+(a011/s);
t21=(b21+(a021/s));
t12=(b12+(a012/s));
t22=b22+(a022/s);

p11=0;
p12=0;
p21=0;
p22=0;

for k=1:2
    Bdelta=((alpha+beta)^2/(alpha^2)/(beta^2))*((((alpha*delta-gamma)*(beta*delta+gamma))/(alpha*beta*(alpha+beta)^2))-((k^2)*(pi^2)/L^2));
    s_Bdelta = sqrt(Bdelta);
        p_1 = -((alpha-beta)*gamma+2*alpha*beta*delta)/((alpha+beta)^2);
        p_2 = (2*(alpha^2)*(beta^2)*s_Bdelta)/((alpha+beta)^2);
    if imag(p_2)==0
        p1 = p_1+p_2;
        p2 = p_1-p_2;
        for p = [p1 p2]
            
            h1=(((alpha-beta)*p+gamma)*L/(2*alpha*beta));
            h2=(T*L^3*p*s_Bdelta);

            a11=-((2*(k^2)*pi^2)/h2);
            a12=(-1)^k*((2*(k^2)*pi^2*exp(-h1))/h2);
            a21=-((2*(k^2)*pi^2*exp(h1))/(h2))*((cos(k*pi))+((h1/(k*pi))*sin(k*pi)));
            a22=(-1)^k*((2*(k^2)*pi^2)/(h2))*((cos(k*pi))+((h1/(k*pi))*sin(k*pi)));
            
            p11 = p11 + (s*a11)/(p*(s-p));
            p12 = p12 + (s*a12)/(p*(s-p));
            p21 = p21 + (s*a21)/(p*(s-p));
            p22 = p22 + (s*a22)/(p*(s-p));
        end
        
    else
        p = p_1 + p_2;
        p_r = real(p);
        p_i = imag(p);
        
        h1=((alpha-beta)*p+gamma)*L/(2*alpha*beta);
        h2=(T*L^3*p*s_Bdelta);
        
        a11 = -((2*(k^2)*pi^2)/h2);
        a12=(-1)^k*((2*(k^2)*pi^2*exp(-h1))/h2)*(cos(k*pi));
        a21=-(-1)^k*((2*(k^2)*pi^2*exp(h1))/h2)*(cos(k*pi));
        a22=((2*(k^2)*pi^2)/(h2));
        
%       real and imaginary parts
        a11_r = real(a11);
        a11_i = imag(a11);
        a12_r = real(a12);
        a12_i = imag(a12);
        a21_r = real(a21);
        a21_i = imag(a21);
        a22_r = real(a22);
        a22_i = imag(a22);
        
         p11 = p11 + 2*s*(s*(a11_r*p_r+a11_i*p_i)-(a11_r*p_r^2-a11_r*p_i^2+2*a11_i*p_r*p_i))/((p_r^2+p_i^2)*(s^2-2*p_r*s+p_r^2+p_i^2));
         
         p12 = p12 + 2*s*(s*(a12_r*p_r+a12_i*p_i)-(a12_r*p_r^2-a12_r*p_i^2+2*a12_i*p_r*p_i))/((p_r^2+p_i^2)*(s^2-2*p_r*s+p_r^2+p_i^2));
         
         p21 = p21 + 2*s*(s*(a21_r*p_r+a21_i*p_i)-(a21_r*p_r^2-a21_r*p_i^2+2*a21_i*p_r*p_i))/((p_r^2+p_i^2)*(s^2-2*p_r*s+p_r^2+p_i^2));
         
         p22 = p22 + 2*s*(s*(a22_r*p_r+a22_i*p_i)-(a22_r*p_r^2-a22_r*p_i^2+2*a22_i*p_r*p_i))/((p_r^2+p_i^2)*(s^2-2*p_r*s+p_r^2+p_i^2));
         
    end

end

P11 = t11 + p11;
P12 = exp(-tu*s)*(t12 + p12);
P21 = exp(-td*s)*(t21 + p21);
P22 = t22 + p22;

xx=[Qm,tu,td];
end


