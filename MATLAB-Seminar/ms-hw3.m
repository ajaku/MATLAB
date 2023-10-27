% MATLAB HW#3 - Andy Jaku
clc
clear
close all

%% Space Between 

%1.1
A = reshape(exp(0:63),8,[]).';
A(3:6,3:6) = A(3:6,3:6)-A(3:6, 3:6);
gmA = geomean(A(A>0));

%1.2
B = sin(linspace(0, 5, 100).*linspace(-5,0, 100).');
[row, column] = find(B > .4995 & B < .5005);
C = diag(diag(B(row,column)));

%1.3
x = -4:.01:4;
y = -4:.01:4
[X,Y] = meshgrid(x,y);
Z_1= exp(-(1-X.*Y).^2);
Z_2= (1/4).*((X.^2.+Y.^2).^(1/2));
S_1=trapz(y,trapz(x,Z_1,2));
S_2=trapz(y,trapz(x,Z_2,2));
S_2 - S_1

%% I Need a Vacation

%2.1
[X,Y] = meshgrid(1:1:256,1:1:256);
A = boolean(((X-99).^(2)+(Y-99).^(2)).^(1/2)<29);
FA = figure;
imshow(A); % Small circle shifted left and slightly up

%2.2
[X,Y] = meshgrid(1:1:256,1:1:256);
B = boolean(((X-62).^(2)+(Y-62).^(2)).^(1/2)<58);
FB = figure;
imshow(B); % Larger circle shifted left and up

%2.3
[X,Y] = meshgrid(1:1:256,1:1:256);
C = boolean((X-(4.*sin(Y./10)))>200);
FC = figure;
imshow(C); % Looks like a wave pattern on the right of graph

%2.4
S = rand(256,256,3);
S(:,:,(1:2)) = 0;
FS = figure;
imshow(S); % Looks like a piece of blue fabric

%2.5
M = A & ~B;
FM = figure;
imshow(M); % Looks like the moon

%2.6
Z = (C.*S) + M;
FZ = figure;
imshow(Z); % Looks like the moon over the ocean, manipulating the tide


%% My Sinc is Broken

%3.1 
%deriv = @(y,x) diff(y);
%antideriv = @(y,x) int(y);
lo = -100;
hi = 100;
N = 1e2;
x = linspace(lo, hi, N);
y = x.^2;
derivVal = deriv(y,x);
integralVal = antideriv(y,x,hi,lo,N);
%Functions are defined at end of script


%3.2
a = [-1 2 1 0 1 -3];
a(a>=0)=1;
a(a<0)=0;
b = circshift(a,1,2);
c = xor(a,b);

%3.3
x_etrema = extrema(y,x)
y_extrema = 0
inflectionpt = inflections(y,x)

%3.4

function [x y] = inflections(a,b)
    minimum = min(min(a))
    [x y] = find(a==minimum)
    maximum = max(max(a))
    [x y] = find(a==maximum)
end
function [x y] = extrema(a,b)
    x = find(deriv(a,b) == 0);
    y = 0;
end
function [val] = deriv(y,x)
    val = diff(y)./diff(x);
end
function [val] = antideriv(y,x,hi,lo,N)
    val = ((hi-lo)/(N-1))*cumsum(y);
end