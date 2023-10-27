% Signals HW#9- Andy Jaku
clc
clear
close all

%% Q-3)
% Part A
E = [3, 1; 2,1];
A_0 = [-0.3, 0.4; -0.4, -0.3];
eigenVals = eig(A_0);
% A value
A = E*A_0*inv(E);
% Confirming same eigenvalues
diff = (abs((eig(A)-eig(A_0))) < 0.0001)

% Part B
s = sym('s');
F = (s*eye(2) -A)^-1;
mat = ilaplace(F);
f = matlabFunction(mat);
f_s = 10;
A_d = f(1/f_s);

% Part C
x_d = zeros(2,101);
x_d0 = [2;1];
for n = 0:100
   x_d(:,n+1) = (A_d^(n))*x_d0;
end

% Part D
x_de = zeros(2,101);
x_d0_de = [2;1];
A_de = eye(2,2) + A*(1/f_s);
for n = 0:100
   x_de(:,n+1) = (A_de^(n))*x_d0_de;
end
x_de;
x_dm = zeros(2,101);
x_d0_dm = [2;1];
A_dm = eye(2,2) + A*(1/f_s) + (A^2)*((1/f_s)^2)*(1/2);
for n = 0:100
   x_dm(:,n+1) = (A_dm^(n))*x_d0_dm;
end

% Part E
hold on
title("State Space Trajectories");
xlabel("x");
ylabel("y");
zlabel("time");
plot3(x_d(1,:), x_d(2,:), (0:0.1:10), "r", displayName = "xd");
plot3(x_de(1,:), x_de(2,:), (0:0.1:10), "b", displayName = "xde");
plot3(x_dm(1,:), x_dm(2,:), (0:0.1:10), "g", displayName = "xdm");
legend();
hold off

% Part F
error_d_de = max(max(abs(x_d - x_de)))
error_d_dm = max(max(abs(x_d - x_dm)))

% Part G

% The midpoint method appears to be more accurate. The error between 
% it and x_d is smaller than x_de, furthermore, the graphs appear to
% overlap substantially more between the x_d and x_dm graphs.

% It is a meaningful amount, especially regarding how similarly related
% the graphs are. x_dm matches x_d quite a bit closer than x_de does.
% The error seems to be substantially larger around curves, which makes
% sense because eulers method can only go in a straight line during a time 
% interval, meanwhile the midpoint method corrects itself once, allowing 
% it to follow curves and bends closer.