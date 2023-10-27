% Signals HW#8- Andy Jaku
clc
clear
close all

%% Q-1) Heavy Tail Distributions
% Part a.
% Generate 1e6 Gaussian, Students (v=5 * factor for var of 1), Cauchy
% (alpha = 0.544)
N = 1e6;
gauss = randn(N,1);
uniform = rand(N,1);
cauchy = (0.544)*tan(pi*uniform);
students = trnd(5, N, 1);
students = 1/sqrt(var(students)) * students;

figure();
histogram(cauchy);
title("Cauchy Distribution");
xlabel("x");
ylabel("y");
xline(-1, "--r");
xline(1, "--b");

figure();
histogram(gauss);
title("Gaussian Distribution");
xlabel("x");
ylabel("y");
xline(-1, "--r");
xline(1, "--b");

figure();
histogram(students);
title("Students Distribution");
xlabel("x");
ylabel("y");
xline(-1, "--r");
xline(1, "--b");

% Find percent of time under 1
gaussTime = sum(abs(gauss(:,1)) < 1)/N;
cauchyTime = sum(abs(cauchy(:,1)) < 1)/N;
studentsTime = sum(abs(students(:,1)) < 1)/N;

% Part b.
changedGauss = reshape(gauss,[100000, 10]);
changedCauchy = reshape(cauchy,[100000, 10]);
changedStudents = reshape(students,[100000, 10]);

gaussMean = mean(changedGauss);
cauchyMean = mean(changedCauchy);
studentsMean = mean(changedStudents);
% Cauchy does not appear to have a mean of 0.

%% Q-2) ARMA and AR Models
% Part A
% Q6
A = [1, -1.6, 0.81]; 
B = [1, 0.4, 0.2]; 
[z,p,k] = tf2zpk(B, A);
figure();
zplane(B, A)
flag = isminphase(B,A); % is minimum-phase

% Part B
% Q1
N = 1e4;
whiteNoise = sqrt(2)*randn(N,1);
x = filter(B,A,whiteNoise);

% Q2
r_x = zeros(7,1);
for m = 0:6
    for n = m+1:N-m
        r_x(m+1) = r_x(m+1) + (x(n+m)*conj(x(n)));
    end
    r_x(m+1) = r_x(m+1)/N;
end

% Q3
figure();
stem(-6:6,r_x(abs(-6:6) +1), "r", "filled");
title("r_x");

% Q4
R = toeplitz(r_x, r_x.');

% Q5
eigenVals = eig(R); % All eigen values are +ve

% Q6
M = 7;
A = zeros(M, N-M+1);
for col = 1:N-M+1
    for row = 1:M
        A(row,col) = x(M-row+col);
    end
end
R_2 = (A*A.')/(N-M+1);
error = norm(R - R_2);

% Part C
% Q1
[s_est, w] = pwelch(x, hamming(512),256,512);
figure();
plot(w, s_est);
title("PSD");
ylabel("Power");
xlabel("Frequency");
vals = [w, s_est];

% Q2
% The maximum value consitently appears just before 0.48, ranging from 0.465 
% to 0.50

% Q3
angle_1 = atan(imag(p(1)/real(p(1))));
angle_2 = atan(imag(p(2)/real(p(2))));
% These values roughly match the estimated value of w_0

% Part 4
[a, varv] = aryule(x,4);
diff = abs(varv -var(whiteNoise)); % Difference consitently appears to be 
% small, suggesting the values are the same.
x_0 = filter(1,a,whiteNoise);
r_x_0 = zeros(7,1);
for m = 0:6
    for n = m+1:N-m
        r_x_0(m+1) = r_x_0(m+1) + (x_0(n+m)*conj(x_0(n)));
    end
    r_x_0(m+1) = r_x_0(m+1)/N;
end
figure();
title("r_x vs rx_0")
hold on
stem(-6:6,r_x(abs(-6:6) +1), "r", "filled");
stem(-6:6,r_x_0(abs(-6:6) +1), "b", "filled");
hold off

figure();
title("x vs x_0")
hold on
stem(x(1:100), "r", "filled")
stem(x_0(1:100), "b", "filled")
hold off
% As was noted in the teams channel, the AR seems to follow the ARMA
% model quite well.