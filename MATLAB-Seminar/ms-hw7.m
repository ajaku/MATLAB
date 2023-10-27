% MATLAB HW#7 - Andy Jaku
clc
clear
close all

% Cheating at Calculus
%% 1 Diff Eq
syms y(x);
cond = y(0) == 2;
ode = diff(y,x,1) == (y^2)*(x^3);
ySol(x) = dsolve(ode, cond);
figure()
fplot(ySol(x))
cond = y(0) == 3;
ySol(x) = dsolve(ode, cond);
figure();
fplot(ySol(x))
cond = y(0) == 100;
ySol(x) = dsolve(ode, cond);
figure();
fplot(ySol(x))

%% 2 Laplace Transform
syms t
h = 5*dirac(t) + heaviside(t) + exp((-t/2))*sin(t)
laplace(h)
ilaplace(h)

%% 3 Vector Fields
syms x y z t a;
F = [((3+a)^2)*(x*z), z*exp(y), (exp(y)-((x^2)*exp(pi*a)))];
r = [(1-2*cos(t))*cos(3*t), (1-2*cos(t))*sin(3*t), sin(t)];
F = subs(F, x, r(1));
F = subs(F, y, r(2));
F = subs(F, z, r(3));
rDiff = diff(r)
integrand = dot(F,rDiff)
work = int(integrand,t,[-pi, pi])
work = subs(work, a, [-3 -2 -1 0 1])
fplot(work)

% It Ain't Magic
%% 1
syms x y
F = 1 - x^2 -y^2
ht = matlabFunction(F)

% 2
% Int method
tic
integral = int(int(F,y,-sqrt(1-x^2),sqrt(1-x^2)),x,-1,1)
toc

tic
x = -1:.1:1;
y = -1:.1:1;
[X,Y] =  meshgrid(x,y);
F_t = 1- X.^2 -Y.^2;
I = trapz(y,trapz(x,F_t,2))
toc