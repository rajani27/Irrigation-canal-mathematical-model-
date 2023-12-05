clear all
clc
% Canal characteristics data
H=3.6;
m=1.5;
n=0.018; 
Sb=0.00025;
L=4846;
B=9.45;
y=1.07;
yu=y;

%Uniform flow case 
yd=y;

%Gate characteristics

G1B=B;
G2B=B;

% Obtain transfer functions and other parameters from UniIDZ
[P11,P12,P21,P22,Qm,xx]=UniIDZ(m,n,Sb,L,B,y);

% Obtain parameters from gategains
Qu=Qm;
Qd=Qm;
[W1,W2,Kw,Ku,Kd,ww]=gategains(Qu,G1B,H,yu,Qd,G2B,yd);