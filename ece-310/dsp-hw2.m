% DSP HW#2- Andy Jaku
clc
clear
close all

%% Q-1)
% Attached to submission 

%% Q-2)
% Specs:
% Creating a bandstop filter
% Pasband below 9MHz and above 12.5MHz
% Passband variation 1.5dB
% Stopband from 9.5MHz to 12MHz
% Stopband attenuation 40db

% For digital filters sampling rate if 40MHz

f_s = 40e6; % Sampling frequency
f_n = f_s/2; % Nyquist frequency

f = linspace(0, f_n, 1001); % Suggested by Logan Po

Wp_dig = [9e6 12.5e6]/(f_n); 
Ws_dig = [9.5e6 12e6]/(f_n);
% Wp(1) = 9 < Ws(1) = 9.5 < Ws(2) = 12 < Wp(2) = 12.5
% As a result, passband should range from 0->Wp(1) and Wp(2) to 1.
Wp_an = [9e6 12.5e6]*(2*pi);
Ws_an = [9.5e6 12e6]*(2*pi);

Rp = 1.5; % passband attenuation
Rs = 40; % stopband attenuation
% if you want analog, rep freq in radians/second and specify 's'.
% As a result, passband should range from 0->Wp(1) and Wp(2) to 1.

%% buttord -> butter (ANALOG)
[n_b_an, Wn_b_an] = buttord(Wp_an, Ws_an, Rp, Rs, 's');
fprintf("\nButterworth Filter Order (Analog): %0.4f\n", 2*n_b_an);
[b_b_an, a_b_an] = butter(n_b_an, Wn_b_an, 'stop', 's');
[z_b_an, p_b_an, ~] = butter(n_b_an, Wn_b_an, 'stop', 's');

% b) Pole-Zero plots 
figure();
zplane(z_b_an, p_b_an);
title('Poles and Zeros: Butterworth Analog');
grid;
h_b_an = freqs(b_b_an, a_b_an, 2*pi*f);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_b_an(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_b_an(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_b_an(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_b_an(f/1e6 == 12))));

figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_b_an)));
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Butterworth Analog');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_b_an))*180/pi);
title('Phase Response: Butterworth Analog');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

% buttord -> butter (DIGITAL)

% a) Filter Order
[n_b_dig, Wn_b_dig] = buttord(Wp_dig, Ws_dig, Rp, Rs);
fprintf("\nButterworth Filter Order (Digital): %0.4f\n", 2*n_b_dig);
[b_b_dig, a_b_dig] = butter(n_b_dig, Wn_b_dig, 'stop');
[z_b_dig, p_b_dig, ~] = butter(n_b_dig, Wn_b_dig, 'stop');

% b) Pole-Zero plots 
figure();
zplane(z_b_dig, p_b_dig);
title('Poles and Zeros: Butterworth Digital');
grid;

% c-d) Magnitude and Phase Response
figure();
h_b_dig = freqz(b_b_dig, a_b_dig, f, f_s);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_b_dig(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_b_dig(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_b_dig(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_b_dig(f/1e6 == 12))));

subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_b_dig)));
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Butterworth Digital');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_b_dig))*180/pi);
title('Phase Response: Butterworth Digital');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

%% cheb1ord -> cheby1(ANALOG)
[n_c1_an, Wp_c1_an] = cheb1ord(Wp_an, Ws_an, Rp, Rs, 's');
fprintf("\nChebyshev I Filter Order (Analog): %0.4f\n", 2*n_c1_an);
[b_c1_an, a_c1_an] = cheby1(n_c1_an, Rp, Wp_c1_an, 'stop', 's');
[z_c1_an, p_c1_an, ~] = cheby1(n_c1_an, Rp, Wp_c1_an, 'stop', 's');

% b) Pole-Zero plots 
figure();
zplane(z_c1_an, p_c1_an);
title('Poles and Zeros: Chebyshev I Analog');
grid;
h_c1_an = freqs(b_c1_an, a_c1_an, 2*pi*f);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_c1_an(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_c1_an(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_c1_an(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_c1_an(f/1e6 == 12))));

figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_c1_an)));
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Chebyshev I Analog');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_c1_an))*180/pi);
title('Phase Response: Chebyshev I Analog');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');


% cheb1ord -> cheby1 (DIGITAL)
% a) Filter Order
[n_c1_dig, Wp_c1_dig] = cheb1ord(Wp_dig, Ws_dig, Rp, Rs);
fprintf("\nChebyshev Type I Filter Order (Digital): %0.4f\n", 2*n_c1_dig);
[b_c1_dig, a_c1_dig] = cheby1(n_c1_dig, Rp, Wp_c1_dig, 'stop');
[z_c1_dig, p_c1_dig, ~] = cheby1(n_c1_dig, Rp, Wp_c1_dig, 'stop');

% b) Pole-Zero plots 
figure();
zplane(z_c1_dig, p_c1_dig);
title('Poles and Zeros: Chebyshev 1 Digital');
grid;

% c-d) Magnitude and Phase Response
%figure();
%title('Chebyshev 1 Bandstop Digital Filter');
%freqz(b_c1_dig, a_c1_dig, [], f_s)
figure();
h_c1_dig = freqz(b_c1_dig, a_c1_dig, f, f_s);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_c1_dig(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_c1_dig(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_c1_dig(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_c1_dig(f/1e6 == 12))));

subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_c1_dig)))
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Chebyshev Type I Digital');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_c1_dig))*180/pi);
title('Phase Response: Chebyshev Type I Digital');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

%% cheb2ord -> cheby2(ANALOG)
[n_c2_an, Wp_c2_an] = cheb2ord(Wp_an, Ws_an, Rp, Rs, 's');
fprintf("\nChebyshev II Filter Order (Analog): %0.4f\n", 2*n_c2_an);
[b_c2_an, a_c2_an] = cheby2(n_c2_an, Rs, Wp_c2_an, 'stop', 's');
[z_c2_an, p_c2_an, ~] = cheby2(n_c2_an, Rs, Wp_c2_an, 'stop', 's');

% b) Pole-Zero plots 
figure();
zplane(z_c2_an, p_c2_an);
title('Poles and Zeros: Chebyshev II Analog');
grid;
h_c2_an = freqs(b_c2_an, a_c2_an, 2*pi*f);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_c2_an(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_c2_an(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_c2_an(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_c2_an(f/1e6 == 12))));

figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_c2_an)));
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Chebyshev II Analog');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_c2_an))*180/pi);
title('Phase Response: Chebyshev II Analog');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');


% cheb2ord -> cheby2 (DIGITAL)
% a) Filter Order
[n_c2_dig, Ws_c2_dig] = cheb2ord(Wp_dig, Ws_dig, Rp, Rs);
fprintf("\nChebyshev Type II Filter Order (Digital): %0.4f\n", 2*n_c2_dig);
[b_c2_dig, a_c2_dig] = cheby2(n_c2_dig, Rs, Ws_c2_dig, 'stop');
[z_c2_dig, p_c2_dig, ~] = cheby1(n_c1_dig, Rs, Wp_c1_dig, 'stop');

% b) Pole-Zero plots 
figure();
zplane(z_c2_dig, p_c2_dig);
title('Poles and Zeros: Chebyshev Type II Digital');
grid;

% c-d) Magnitude and Phase Response
h_c2_dig = freqz(b_c2_dig, a_c2_dig, f, f_s);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_c2_dig(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_c2_dig(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_c2_dig(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_c2_dig(f/1e6 == 12))));

figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_c2_dig)))
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Chebyshev Type II Digital');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_c2_dig))*180/pi)
title('Phase Response: Chebyshev Type II Digital');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

%% ellipord -> ellip (ANALOG)
[n_e_an, Wp_e_an] = ellipord(Wp_an, Ws_an, Rp, Rs, 's');
fprintf("\nElliptical Filter Order (Analog): %0.4f\n", 2*n_e_an);
[b_e_an, a_e_an] = ellip(n_e_an, Rp, Rs, Wp_e_an, 'stop', 's');
[z_e_an, p_e_an, ~] = ellip(n_e_an, Rp, Rs, Wp_e_an, 'stop', 's');


% b) Pole-Zero plots 
figure();
zplane(z_e_an, p_e_an);
title('Poles and Zeros: Elliptical Analog');
grid;
h_e_an = freqs(b_e_an, a_e_an, 2*pi*f);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_e_an(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_e_an(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_e_an(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_e_an(f/1e6 == 12))));

figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_e_an)));
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Elliptical Analog');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_e_an))*180/pi);
title('Phase Response: Elliptical II Analog');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

% ellipord -> ellip (DIGITAL)
% a) Filter Order
[n_e_dig, Wp_e_dig] = ellipord(Wp_dig, Ws_dig, Rp, Rs);
fprintf("\nElliptical Filter Order (Digital): %0.4f\n", 2*n_e_dig);
[b_e_dig, a_e_dig] = ellip(n_e_dig, Rp, Rs, Wp_e_dig, 'stop');
[z_e_dig, p_e_dig, ~] = ellip(n_e_dig, Rp, Rs, Wp_e_dig, 'stop');

% b) Pole-Zero plots 
figure();
zplane(z_e_dig, p_e_dig);
title('Poles and Zeros: Elliptical Digital');
grid;

% c-d) Magnitude and Phase Response
figure();
h_e_dig = freqz(b_e_dig, a_e_dig, f, f_s);
fprintf("Peak Passband @ 9MHz Edge: %0.4f\n", 20*log10(abs(h_e_dig(f/1e6 == 9))));
fprintf("Peak Passband @ 12.5MHz Edge: %0.4f\n", 20*log10(abs(h_e_dig(f/1e6 == 12.5))));
fprintf("Stopband Edge @ 9.5MHz: %0.4f\n", 20*log10(abs(h_e_dig(f/1e6 == 9.5))));
fprintf("Stopband Edge @ 12MHz: %0.4f\n", 20*log10(abs(h_e_dig(f/1e6 == 12))));
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_e_dig)))
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Elliptical Digital');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_e_dig))*180/pi)
title('Phase Response: Elliptical Digital');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

%% Q-3)
% Specs:
% Creating a bandstop filter
% Pasband below 9MHz and above 12.5MHz
% Passband variation 1.5dB
% Stopband from 9.5MHz to 12MHz
% Stopband attenuation 40db

% For digital filters sampling rate if 40MHz

f_s = 40e6; % Sampling frequency
f_n = f_s/2; % Nyquist frequency

f = linspace(0, f_n, 1001); % Suggested by Logan Po

Wp_dig = [9e6 12.5e6]/(f_n); 
Ws_dig = [9.5e6 12e6]/(f_n);
% Wp(1) = 9 < Ws(1) = 9.5 < Ws(2) = 12 < Wp(2) = 12.5
% Passband Gain = 20log10((1+delta_passband)/(1-delta_passband))
% Stopband Gain = -20log10(delta_stopband)
% Solve for both (known that Passband Gain = 1.5dB, Stopband = 40dB)
syms x
eqn = 20*log10((1+x)/(1-x)) == 1.5;
delta_passband = eval(solve(eqn, x));
fprintf("\nDelta Passband: %0.4f\n", delta_passband);
syms x
eqn = -20*log10(x) == 40;
delta_stopband = eval(solve(eqn, x));
fprintf("Delta Stopband: %0.4f\n", delta_stopband);


% Kaiser 
[n_kai, Wn_kai, beta_kai, ftype_kai] = kaiserord([9e6 9.5e6 12e6 12.5e6], [1 0 1], [delta_passband delta_stopband delta_passband], f_s);
fprintf("\nKaiser Filter Order (Digital): %0.4f\n", 2*n_kai);
b_kai = fir1(n_kai,Wn_kai,ftype_kai,kaiser(n_kai+1,beta_kai),'noscale');
% b) Pole-Zero plots 
figure();
zplane(b_kai);
title('Poles and Zeros: Kaiser');
grid;

h_kai = freqz(b_kai, 1, f, f_s);

figure();
impz(b_kai,1)
title("Coefficients of Kaiser")


figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_kai)))
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Kaiser');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_kai))*180/pi)
title('Phase Response: Kaiser');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

% Equiripple
[n_eq, f_eq, ao_eq, w_eq] = firpmord([9e6 9.5e6 12e6 12.5e6], [1 0 1], [delta_passband delta_stopband delta_passband], f_s);
fprintf("\nEquiripple Filter Order (Digital): %0.4f\n", 2*n_eq);
b_eq = firpm(n_eq, f_eq, ao_eq, w_eq);
% b) Pole-Zero plots 
figure();
zplane(b_eq, 1);
title('Poles and Zeros: Equiripple');
grid;

h_eq = freqz(b_eq, 1, f, f_s); % suggested by Ridwan

figure();
impz(b_eq,1);
title("Coefficients of Equiripple")


figure();
subplot(2, 1, 1);
plot(f/1e6, 20*log10(abs(h_eq)))
ax = gca;
ax.YLim = [-50 2];
title('Magnitude Response: Equiripple');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
subplot(2, 1, 2);
plot(f/1e6, unwrap(angle(h_eq))*180/pi)
title('Phase Response: Equiripple');
xlabel('Frequency (MHz)');
ylabel('Phase (Degrees)');

% c)
w_eq;
% I think I capped here
fprintf("\nFirpmord Wp/Ws: %0.4f\n", w_eq(1)/w_eq(2));
fprintf("Calculated Ds/Dp: %0.4f\n", delta_stopband/delta_passband);
% It exists

% d)
stop_bound_1 = find(f==9.5e6);
stop_bound_2 = find(f==12e6);
pass_bound_1 = find(f==9e6);
pass_bound_2 = find(f==12.5e6);


fprintf("\nMaximum Passband Kaiser (dB): %0.4f\n", max((20*log10(abs(h_kai(1:pass_bound_1))))) - min((20*log10(abs(h_kai(1:pass_bound_1))))));
fprintf("Maximum Stopband Kaiser (dB): %0.4f\n", min(abs(20*log10(abs(h_kai(stop_bound_1:stop_bound_2))))));
fprintf("Maximum Passband Equiripple (dB): %0.4f\n", max((20*log10(abs(h_eq(1:pass_bound_1))))) - min((20*log10(abs(h_eq(1:pass_bound_1))))));
fprintf("Maximum Stopband Equiripple (dB): %0.4f\n", min(abs(20*log10(abs(h_eq(stop_bound_1:stop_bound_2))))));


% e) Neither fail to meet expectations, but the Kaiser filter is
% overdesigned as shown by its maximum passband being far lower than the
% tolerance specified. To fix that problem you could adjust the tolerances
% for kaiser and make it more lenient.