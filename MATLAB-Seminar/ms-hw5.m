% MATLAB HW#5 - Andy Jaku
clc
clear
close all

%% 1
% 1.1 
num = [3/7 2/3 1/2]; % 1/2 + 2/3z + 3/7z^2
den = [1/2 0 1/3 2]; % 2 + 1/3z +1/2z^3

% 1.2 
figure; 
zplane(num, den);


% 1.3
[Z, P, K] = tf2zp(num,den);

% 1.4
figure;
ir = impz(den, num, 66);
stem(ir);
title('Impulse Response')
xlabel('Impulse (0-66)')
ylabel('Value of Impulse')

% 1.5
sinVals = sin((1:55)/10);
filterVals = filter(num, den, sinVals);
convVals = conv(ir, sinVals);

figure;
subplot(2,1,1);
stem(filterVals);
title('Filtered Sin Values')
xlabel('Sample Points (1-55)')
ylabel('Output of Filter')

subplot(2,1,2);
stem(convVals);
title('Convoloved Sin Values')
xlabel('Sample Points (1-55)')
ylabel('Output of Convolution')

%% Pole-y Moley
% 2.1
load handel.mat

% 2.2
Z = [exp(2.06*j) exp(-2.06*j) exp(1.43*j) exp(-1.43*j)]';
k = .1;
P = [.74*exp(.76*j) .74*exp(-.76*j) .96*exp(1.24*j) .96*exp(-1.24*j)];
[den2, num2] = zp2tf(Z,P,k);

% 2.3 
% Disclaimer: shamelessly pirating code from the lecture note
n = 1024;     % number of samples in frequency response
fs = 20000;   % sampling frequency of input signal
[H, f] = freqz(num2, den2, n, fs);

Hdb = 20*log10(abs(H));
Hph = rad2deg(unwrap(angle(H)));

figure;
subplot(2, 1, 1);
plot(f, Hdb);
xlabel("Frequency (Hz)");
ylabel("|H| (dB)")
xlim([0 fs/2]);
title("Magnitude Response");

subplot(2, 1, 2);
plot(f, Hph);
xlabel("Frequency (Hz)");
ylabel("Phase (deg)");
xlim([0 fs/2]);
title("Phase Response");

% 2.4
soundsc(Hdb)
soundsc(Hph)