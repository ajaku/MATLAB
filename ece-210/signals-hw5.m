% Signals Homework #5 MATLAB - Andy Jaku
clc
clear
close all

%% 1-a
n = 6;
f_lo = 10*10^3;
f_hi = 12*10^3;
r_p = 2;
f_s = 40*10^3;

f_n = f_s/2;
f_crit = [f_lo/f_n, f_hi/f_n];
n0 = n/2;

[z,p,k] = cheby1(n0, r_p, f_crit);
figure;
zplane(z,p);
[b,a] = zp2tf(z,p,k);

% 1-b
f = linspace(0, f_n, 10^4);
[H,~] = freqz(b,a, 1e4);

% 1-c
Hdb = 20*log10(abs(H));
Hph = rad2deg(unwrap(angle(H)));

figure();
subplot(2,1,1);
plot(f/1e3, Hdb);
title("1c) - Magnitude Response");
xlabel("Frequency (kHz)"+newline+"    ");
ylabel("Magnitude (db)");
ylim([-40 2]);
xlim([0 f_n/1e3]);
grid on

subplot(2,1,2)
plot(f/1e3, Hph);
title("1c) - Phase Response")
xlabel("Frequency (kHz)"+newline+"    ")
ylabel("Phase (deg)")
xlim([0 f_n/1e3])
grid on

% 1-d
% Find indicies of values >= -30db
indicies = find(Hdb >= -30);
% Indicies of freq right before beginning and end 
% represent changing points from -29.9999 -> -30.0000
f_low = f(indicies(1)-1);
f_high = f(indicies(end)+1);

%% 4-a
K_o = 4-2*sqrt(2);
K_max = 4;

K_1 = K_o;
K_2 = 0.5*K_o + 0.5*K_max;
K_3 = 0.2*K_o + 0.8*K_max;

R = 10^3;
C = 10*10^-9;
N = 10^4;

[Hdb_1, Hph_1, f_1] = plotQ4a(K_1, R, C, "K_1", -35, -7);
[Hdb_2, Hph_2, f_2] = plotQ4a(K_2, R, C, "K_2", -28, 6);
[Hdb_3, Hph_3, f_3] = plotQ4a(K_3, R, C, "K_3", -26, 17);

% 4-b
generateBodePlots(K_1, R, C, "B_1", -35, -7);
generateBodePlots(K_2, R, C, "B_2", -28, 6);
generateBodePlots(K_3, R, C, "B_3", -26, 17);

% 4-c
[W_1, Q_1] = compWnQ(K_1, R, C);
[W_2, Q_2] = compWnQ(K_2, R, C);
[W_3, Q_3] = compWnQ(K_3, R, C);

% 4-d
[f_peak_1, f_low_1, f_high_1, error_1] = findFrq(Hdb_1, f_1);
[f_peak_2, f_low_2, f_high_2, error_2] = findFrq(Hdb_2, f_2);
[f_peak_3, f_low_3, f_high_3, error_3] = findFrq(Hdb_3, f_3);

% 4-e
[check_1] = checkPeak(f_peak_1, W_1, error_1);
[check_2] = checkPeak(f_peak_2, W_2, error_2);
[check_3] = checkPeak(f_peak_3, W_3, error_3);

% 4-f
[bw_1, f_ctr_1] = findBandwidth(f_low_1, f_high_1);
[bw_2, f_ctr_2] = findBandwidth(f_low_2, f_high_2);
[bw_3, f_ctr_3] = findBandwidth(f_low_3, f_high_3);

% 4-g
[ch_1] = check2ndOrder(f_peak_1, f_ctr_1, Q_1, f_low_1, f_high_1, error_1);
[ch_2] = check2ndOrder(f_peak_2, f_ctr_2, Q_2, f_low_2, f_high_2, error_2);
[ch_3] = check2ndOrder(f_peak_3, f_ctr_3, Q_3, f_low_3, f_high_3, error_3);

% Functions
function [Hdb_a, Hph_a, f_a] = plotQ4a(K, R, C, val, min, max)
    n_a = 10^4;
    f_a = linspace(0,1e6, n_a);
    b_a = [(K/(R*C)) 0];
    a_a = [1 ((4-K)/(R*C)) (2/(R^2*C^2))];
    radianFrq_a = 2*pi*f_a;
    [H_a,~] = freqs(b_a,a_a,radianFrq_a);

    Hdb_a = 20*log10(abs(H_a));
    
    figure;
    subplot(2,1,1);
    plot(f_a/1e3, Hdb_a);
    title("4a) - Magnitude Response "  + val);
    xlabel("Frequency (kHz)"+newline+"    ");
    ylabel("Magnitude (db)");
    ylim([min max]);
    xlim([0 1e6/1e3]);
    grid on
    
    Hph_a = rad2deg(unwrap(angle(H_a)));
    subplot(2,1,2);
    plot(f_a/1e3,Hph_a);
    title("$a) - Phase Response " + val);
    xlabel("Frequency (kHz)"+newline+"    ");
    ylabel("Phase (deg)");
    xlim([0 1e6/1e3]);
    grid on;
end

function generateBodePlots(K, R, C, val, min, max)
    n_b = 10^4;
    f_b = logspace(3,6, n_b);
    b_b = [(K/(R*C)) 0];
    a_b = [1 ((4-K)/(R*C)) (2/(R^2*C^2))];
    radianFrq_b = 2*pi*f_b;
    [H_b,~] = freqs(b_b,a_b,radianFrq_b);

    Hdb_b = 20*log10(abs(H_b));
    
    figure;
    subplot(2,1,1);
    semilogx(f_b/1e3, Hdb_b);
    title("4b) - Magnitude Response "  + val);
    xlabel("Frequency (kHz)"+newline+"    ");
    ylabel("Magnitude (db)");
    ylim([min max]);
    xlim([0 1e6/1e3]);
    grid on;
    
    Hph_b = rad2deg(unwrap(angle(H_b)));

    subplot(2,1,2);
    semilogx(f_b/1e3, Hph_b);
    title("4b) - Phase Response " + val);
    xlabel("Frequency (kHz) "+newline+"    ");
    ylabel("Phase (deg)");
    xlim([0 1e6/1e3]);
    grid on;
end

function [W_n, Q] = compWnQ(K, R, C)
    W_n = sqrt(2)/(R*C);
    Q = sqrt(2)/(4-K);
end

function [f_peak, f_low, f_high, error] = findFrq(Hdb, f)
    peakDb = max(Hdb);
    indicies = find(peakDb == Hdb);
    f_peak = f(indicies);

    idxM3 = find(Hdb > peakDb - 3);
    f_low = f(idxM3(1) - 1);
    f_high = f(idxM3(end) + 1);

    error = f(end) - f(end-1);
end

function [val] = checkPeak(f_peak, W_n, error)
    % Define error as spacing between points on linspace
    % Will say value is valid if < error margin
    val = (2*pi*f_peak - W_n);
    if abs(val) < error
        val = 1;
    else
        val = 0;
    end
end

function [bandwidth, f_ctr] = findBandwidth(f_lo, f_hi)
    bandwidth = f_hi - f_lo;
    f_ctr = sqrt(f_lo*f_hi);
end

function [val] = check2ndOrder(f_peak, f_ctr, Q, f_lo, f_hi, error)
    diff1 = f_peak - f_ctr;
    diff2 = Q - (f_ctr/(f_hi-f_lo));
    if (abs(diff1) < error) && (abs(diff2) < error)
        val = 1;
    else
        val = 0;
    end
end