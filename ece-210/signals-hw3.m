% Signals Homework #3 MATLAB - Andy Jaku

%% 3-a
% test set to see if functions work
x = [1 2 3 4 5 6 7]
M = 4
y = upsample(x,M);
z = upsamp(x,M);
f = downsample(x, M);
g = downsamp(x, M);

%% 3-b
x = 1:7
xup = upsamp(x,2);
xdown = downsamp(x,2);
subplot(3,1,1)
stem((0:length(x)-1),x, "filled")
xlabel("n"+newline+"    ")
ylabel("magnitude")
title("X Standard")

subplot(3,1,2)
stem((0:length(xup)-1),xup,"filled")
title("X upsampled")
xlabel("n"+newline+"    ")
ylabel("magnitude")

subplot(3,1,3)
stem((0:length(xdown)-1),xdown,"filled")
title("X downsampled")
xlabel("n"+newline+"    ")
ylabel("magnitude")

%% 3-c
h = [3, 2, 1, 4];
g1 = upsamp(h, 2);
g2 = downsamp(conv(g1,x),2);
g3 = conv(h,downsamp(x,2));

%% 3-d
difference1 = max(abs(g3-g2));

%% 3-e
noble1 = upsamp(conv(h,x),2);
noble2 = conv(g1,upsamp(x,2));
difference2 = max(abs(noble1 - noble2));

function [a] = upsamp(x,M)
    a = zeros(1,(length(x)-1)*(M)+1);
    a(1:M:end) = x;
end

function [a] = downsamp(x,M)
    a = x(1:M:end);
end