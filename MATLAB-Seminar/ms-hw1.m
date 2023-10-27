% Start
clc
clear
close all
% Scalers
a = abs(sin((pi/3)) + j/(sec((-5*pi)/3)))
l = 8^(1/3)
    syms k
u = (2*symsum(k/factorial(6), k, 1, 80))^(1/2)
m = (imag(floor(log(66^(3.5j)))))^2
% Mother
A={a;l;u;m}
F=[a l; u m]
T=F.'
B=(T*F)^(-1)
C = [T F; F T]
% Cruelty
M=mean(B, 'all')
S=mean(C).'
% Odd Types
D=T+F
E=T+1
F=C+A 
% I think this value makes sense, but I am not brushed up on my matrix
% addition skills

% Not what it means
% Index is 3
INDEX=3
G=[linspace(0,10,INDEX)];
H=G.^2;
I=sum(H);
J=I/INDEX

% Index is 5
INDEX=5
G=[linspace(0,10,INDEX)];
H=G.^2;
I=sum(H);
J=I/INDEX

% Index is 10
INDEX=10
G=[linspace(0,10,INDEX)];
H=G.^2;
I=sum(H);
J=I/INDEX

% Index is 300
INDEX=300
G=[linspace(0,10,INDEX)];
H=G.^2;
I=sum(H);
J=I/INDEX

% Index is 1000000
INDEX=1000000
G=[linspace(0,10,INDEX)];
H=G.^2;
I=sum(H);
J=I/INDEX
%The value appears to be converging towards 33.3333333