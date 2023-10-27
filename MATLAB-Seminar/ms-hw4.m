% MATLAB HW#4 - Andy Jaku
clc
clear
close all

%% Shooting the Schmidt 

%1.1
A = [0 2 7j 8; 9 1 -4 -6j; 6 0 2 -1];
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
V = gramschmidt(A);

%1.2
% identity = is_orthonormal(V);
% Unfortunately, this does not seem to work for 
% complex valued inputs

%1.3
B = [3; 4; 5];
C = ortho_proj(B,V);

%% De-Gauss, Re-Gauss

%2.1
x = linspace(0,2*pi,1000);
sinVals = sin(x);
sigma = 1;

mew = 0;
gaus1 = (1/((2*pi*(sigma.^(2)).^(1/2)))*exp((-(x-mew).^2)/(2*(sigma.^2))));

mew = (pi/2);
gaus2 = (1/((2*pi*(sigma.^(2)).^(1/2)))*exp((-(x-mew).^2)/(2*(sigma.^2))));

mew = (pi);
gaus3 = (1/((2*pi*(sigma.^(2)).^(1/2)))*exp((-(x-mew).^2)/(2*(sigma.^2))));

mew = (3*pi/2);
gaus4 = (1/((2*pi*(sigma.^(2)).^(1/2)))*exp((-(x-mew).^2)/(2*(sigma.^2))));

mew = (2*pi);
gaus5 = (1/((2*pi*(sigma.^(2)).^(1/2)))*exp((-(x-mew).^2)/(2*(sigma.^2))));

V = [gaus1; gaus2; gaus3; gaus4; gaus5]'

%2.2
set = gramschmidt(V);
ortho = is_orthonormal(set);
sinVals
sinProj = ortho_proj(sinVals, set);
subplot(2,1,1);
plot(x, sinVals)
title('sin(x)')
xlabel('0 - 2pi')
ylabel('Value of function')

subplot(2,1,2);
plot(x, sinProj)
title('Projected sin(x)')
xlabel('0 - 2pi')
ylabel('Value of function')


function U = gramschmidt(V)
    [n, k] = size(V);
    U = zeros(n,k);
    U(:,1) = V(:,1) / norm(V(:,1)); % Special case for first column (W1 = V1)
    for i = 2:k
        U(:,i) = V(:,i);
        for j = 1:i-1
            U(:,i) = U(:,i) - ( dot(U(:,i),U(:,j)) )*( U(:,j) )/( norm(U(:,j))^(2) );
            % U(:,i) = V (n)
            % U(:,j) = W (opposite direction)
        end
        U(:,i) = U(:,i) / norm(U(:,i));
        % Normalize for phi values
    end
end
% Credit: I altered this algorithim to work with how we learned in class.
% however, the original algorithm idea was based on this:
% https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process

function U = is_orthonormal(V)
    A = V'*V; % will give identity matrix if orthonormal
    I = eye(rank(A));
    B = A-I;
    items = find(B<0.0000001);
    if numel(items) == numel(B)
        U = 1;
    else
        U = 0;
    end
end

function U = ortho_proj(C,V)
    [n,k] = size(V);
    U = zeros(n,1);
    for j = 1:k
        U = U + dot(C,V(:,j))/dot(V(:,j),V(:,j))*V(:,j); %Ortho proj formula looped over columns
    end
end