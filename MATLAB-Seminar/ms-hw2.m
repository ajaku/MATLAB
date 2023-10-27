% MATLAB HW#2 - Andy Jaku
clc
clear
close all

%vNAME = vector named NAME
%mNAME = matrix named NAME

%% Space 

vA = exp(0:63);
mA = reshape(vA, 8, 8, []);
mA = transpose(mA);
mB = mA(3:6,3:6); % Starting idx + 1 : end idx
mA(3:6, 3:6) = mA(3:6,3:6) - mB;
gmB = (prod(mB, "all"))^(1/numel(mB)); % Geometric mean of B
mA(:,1) = flipud(mA(:, 1)); 
mA([4,5],:) = [];


%% Speed (Longest Run)

a = linspace(0, 1000*pi, 10000);
b = linspace(0, 1000*pi, 10000);
tic
 for a =1:10000
     for b = 1:10000
         mC(a,b)=sin(abs(a*(1000*pi)+j*b*(1000*pi)))/(abs(a*(1000*pi)+j*b*(1000*pi)));
     end
 end
toc

% Took 371.573197 seconds (6.183 minutes)

%% Speed (Second Longest)

tic
mD=zeros(10000);
for a =1:10000
    for b = 1:10000
        mD(a,b)=sin(abs(a*(1000*pi)+j*b*(1000*pi)))/(abs(a*(1000*pi)+j*b*(1000*pi)));
    end
end
toc

% Took 34.243477 seconds

%% Speed (Second Shortest)

tic
a = linspace(0, 1000*pi, 10000);
b = linspace(0, 1000*pi, 10000);
[X,Y]=meshgrid(a,b);
F=(sin(abs(X+j*Y))./abs(X+j*Y));
toc

%Took 0.898946 seconds

%% Speed (Shortest)

tic
a = linspace(0, 1000*pi, 10000);
b = linspace(0, 1000*pi, 10000);
b = transpose(b);
mE = (sin(abs(a+j*b))./abs(a+j*b));
toc

%Took 0.719931
 
surf(mCCC(1:50,1:50));

%% The Long Sssss

x1 = linspace(0, 6.66, 100);
x2 = linspace(0, 6.66, 10000);
x1 = x1.*exp(-x1.^2);
x2 = x2.*exp(-x2.^2);

dx1 = 6.66/(100-1)
dx2 = 6.66/(10000-1)

y1 = (2/((pi)^(1/2)))*cumsum(x1).*dx1;
y2 = (2/((pi)^(1/2)))*cumsum(x2).*dx2;

y3 = (2/((pi)^(1/2)))*cumtrapz(x1).*dx1;
y4 = (2/((pi)^(1/2)))*cumtrapz(x2).*dx2;


mse_y1 = mean((y1 - erf(x1)) .^ 2);
mse_y2 = mean((y2 - erf(x2)) .^ 2);
mse_y3 = mean((y3 - erf(x1)) .^ 2);
mse_y4 = mean((y4 - erf(x2)) .^ 2);

plot(y4);