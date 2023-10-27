% DSP HW#1- Andy Jaku
clc
clear
close all

%% Q-1)
% 1a) Nyquist bandwidth, f_s/2 => 20MHz/2 = 10MHz
% 1b) 6MHz, 14MHz, 26MHz, 34MHz
% 1c)
%   Angular Radian Freq: w=2*pi*f = 12*piMHz = 12*piM(rad/s)
%   Digital Radian Freq: w=2*pi*f/fs = 12*pi/20 = 3/5*pi (rad)
%   Digital Freq Normalized to Sampling Rate: f/fs = 6/20 = 3/10
%   Digital Freq Normalized to Nyquist Bandwidth: f/B_nyq = 6/10 = 3/5
% 1d) 6MHz, 44MHz, 56MHz, 94MHz

%% Q-2)
% 2a)
f = 6*10^6;
f_s = 20*10^6;
N_samp = 500;
T = 1/f_s; % Calculate the period needed for 512 samples
t_temp = (0:N_samp-1)*T;
t(1, 1:512) = 0;
t(1, (1:500)) = t_temp;

S = sin(2*pi*f*t_temp);
N = 512;
X = fft(S, 512);

[B, I] = maxk(X, 4); %Maxes occur at 154 and 358


X_shift = fftshift(X);

k = (0:N-1);
f_freq = k*f_s/N; % analog freq (MHz)
w = 2*pi*k/N;     % digital rad freq (rads)

figure(1);
plot((0:N-1), abs(X))
title("2b) FFT Magnitude vs N-th Sample")
xlabel("N Value")
ylabel("|fft(X)|")

plot(f_freq, abs(X))
title("2b) FFT Magnitude vs Frequency")
xlabel("f (Hz)")
ylabel("|fft(X)|")

plot(w, abs(X))
title("2b) FFT Magnitude vs Digital Rad Freq")
xlabel("radians (rad)")
ylabel("|fft(X)|")

% Graph only the positive segment
k_pos = 0:N/2;
X_pos = X(1:N/2+1);
f_freq = k_pos*f_s/N; % analog freq (MHz)
plot(f_freq, abs(X_pos))
title("2b) FFT Magnitude vs Positive Frequency")
xlabel("f (Hz)")
ylabel("|fft(X)|")

% 2c)

% Compute DFT using rectangular window
rect_win = rectwin(N_samp)';
rect_energy_scale = norm(rect_win);
rect_win = rect_win/rect_energy_scale;

X_rect = fft(S.*rect_win, 512);

k_o = 154;

% Compute DFT using chebyshev window
cheb_win = chebwin(N_samp, 30)';
cheb_energy_scale = norm(cheb_win);
cheb_win = cheb_win/cheb_energy_scale;


X_cheb = fft(S.*cheb_win, 512);

figure(3)
plot((k_o-10):(k_o+10), abs(X_rect(k_o-10:k_o+10)), 'Color', 'blue')
title("2c) Rectangular and Chebyshev Magniutde vs N")
xlabel("N Value")
ylabel("Magnitude of DFT");
hold on
plot((k_o-10):(k_o+10), abs(X_cheb(k_o-10:k_o+10)), 'Color', 'red')
hold off
legend('Rectangular','Chebyshev')



% 2d)
snr_level = 20;
S_noise = 10^((-snr_level)/20)*randn(1, 500) + S;

figure(4);
plot((0:500-1), S_noise)
title("2d) Noisey Sine Wave vs N-th Sample")
xlabel("N Value")
ylabel("Sine Func")

X_noise_rect = fft(S_noise.*rect_win, 512);
X_noise_rect_shift = fftshift(X_noise_rect);

X_noise_cheb = fft(S_noise.*cheb_win, 512);
X_noise_cheb_shift = fftshift(X_noise_cheb);

k = (-N/2:N/2-1);
f_freq = k.*f_s/N; % analog freq (MHz) NOT SURE ABOUT ELEMENT WISE MULT
w = 2*pi*k/N;     % digital rad freq (rads)

figure(5);
plot(f_freq, 20*log10(abs(X_noise_rect_shift)))
title("2d) Decibel Magnitude of Noisey Rectangular DFT vs. Frequency")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")

figure(6);
plot(f_freq, 20*log10(abs(X_noise_cheb_shift)))
title("2d) Decibel Magnitude of Noisy Chebyshev DFT vs. Frequency")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")


[B, I] = maxk(20*log10(abs(X_noise_rect_shift)), 2); % Peaks now occur at N= 103, and N = 411
[B, I] = maxk(20*log10(abs(X_noise_cheb_shift)), 6); % Peaks now occur at N= 103, and N = 411
% However, applying new shifted K index, these values should be subtracted
% by N/2 to shift them to match the indexing of the shifted FFT.
% Thus, the N values are 103-256 = -153 and 411-256 = 155, which match the
% expected values from the graph of N/2:N/2-1
bin_spacing = f_s/512; %0.04*10^6;
max_of_shift = 103;
shift_scale = 256;
corrected_val = shift_scale - max_of_shift;
f_o = corrected_val*bin_spacing;

figure(7)
plot(((f_o-10*bin_spacing):bin_spacing:(f_o+10*bin_spacing)), 20*log10(abs(X_noise_rect_shift(max_of_shift-10:max_of_shift+10))), 'Color', 'blue')
title("2d) Noisey Rectangular and Chebyshev Magniutde vs Frequency")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)");
hold on
plot(((f_o-10*bin_spacing):bin_spacing:(f_o+10*bin_spacing)), 20*log10(abs(X_noise_cheb_shift(max_of_shift-10:max_of_shift+10))), 'Color', 'red')
hold off
legend('Rectangular','Chebyshev')

%% Q-3)
% a)
f_s = 44.1e3;     % 44.1kHz
bin_spacing = 2;  % Hz
N_0 = f_s/bin_spacing; 
beats = 100;
time_per_beat = 60;
N = N_0*(101/2);
samples_per_beat = time_per_beat/beats*f_s;

f_G4 = 392;
f_A4 = 440;
f_D5 = 587.33;

signal_vec = zeros(1, N);
snr_music = 40;
w_noise = 10^((-snr_music)/20)*randn(1, N);

% notes could have changed after 0.6 seconds thus portion of notes must be
% updated to reflect changed value: Credit @Kengbo Lu for explanation 
for i = 1:samples_per_beat:N
    rv = [1, 1, 1];
    rv(randsample(1:3, 1)) = 0;
    signal_vec(i)=sin((2*pi*f_G4)*(i)/f_s)*rv(1) +sin((2*pi*f_A4)*(i)/f_s)*rv(2) + sin((2*pi*f_D5)*(i)/f_s)*rv(3) + w_noise(i);
    % attempting to account for beat overlap due to tempo
    for j = i:N-i-samples_per_beat
        signal_vec(j)=sin((2*pi*f_G4)*(j)/f_s)*rv(1) +sin((2*pi*f_A4)*(j)/f_s)*rv(2) + sin((2*pi*f_D5)*(j)/f_s)*rv(3) + w_noise(j);
    end
end


% 3b)
next_power = 2^(nextpow2(N_0));
[pxx,freq]=pwelch(signal_vec, hamming(N_0), next_power/2, next_power, f_s);

% 3c)
figure();
plot(freq,pxx);
xlabel('3c) Frequency (Hz)');
ylabel('PSD Estimate (Watts/Hz)');
title('Periodogram (Hamming window)')

[B, I] = maxk(pxx, 10); %Maxes occur from 290 - 440


figure();
plot(freq(250:500),pxx(250:500));
xlabel('3c) Frequency (Hz)');
ylabel('PSD Estimate (Watts/Hz)');
title('Periodogram (300-600Hz)')

% 3d)
figure;
spectrogram(signal_vec, hamming(N_0), next_power/2, next_power, f_s)
title('3d) Spectrogram (Full Range)')


figure;
spectrogram(signal_vec, hamming(N_0), next_power/2, next_power, f_s)
title('3d) Spectrogram (300-600Hz)')
xlim([0.3, 0.6]);