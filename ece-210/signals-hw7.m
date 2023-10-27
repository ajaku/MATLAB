% Signals HW#7- Andy Jaku
clc
clear
close all

% The Signal
%% Part 1 - Signal Generation
M = 100;
K = 20;
L = 3;
N = 200;
PdB = [0, -2, -4];
PndB = 10;
[A, S] = generateAVec(M, N, K, PdB, PndB);
R = (1/N)*(A*A');

% Part 2 - Analysis
[U, sval, V] = svd(A);
sval = diag(sval);
[eigvec,eigval0]= eig(R);
[eigval,idx]= sort(diag(eigval0),'descend');
eigvec = eigvec(:,idx);

figure();
stem(sval)
title("Singular Values of A (Part II)")

figure();
stem(eigval)
title("Eigen Values of R (Part II)")

ratio_singVal = sval(3)/sval(4);
ratio_eigVal = eigval(3)/eigval(4);
U_L = U(:, 1:L);
P_S = U_L * U_L.';
P_N = eye(length(P_S), length(P_S)) - P_S;
invR = inv(R);

S_rand = generateRandSVec(100, 20, 20);

S_MVDR_Known = genSpectrums(S, invR);
[MVDR_Known_Stats] = [max(S_MVDR_Known), mean(S_MVDR_Known), median(S_MVDR_Known)]

S_MVDR_Rand = genSpectrums(S_rand, invR);
[MVDR_Rand_Stats] = [max(S_MVDR_Rand), mean(S_MVDR_Rand), median(S_MVDR_Rand)]

S_Music_Known = genSpectrums(S, P_N);
[Music_Known_Stats] = [max(S_Music_Known), mean(S_Music_Known), median(S_Music_Known)]

S_Music_Rand = genSpectrums(S_rand, P_N);
[Music_Rand_Stats] = [max(S_Music_Rand), mean(S_Music_Rand), median(S_Music_Rand)]

% Music appears to be better than MVDR by a factor of ~30x
% Furthermore, both MVDR and Music are substantially larger than their
% repsective random values, roughly by a factor of ~6.5x (MVDR) and ~13x (MUSIC) 

%% Part 3 - Try Again
M = 100;
K = 20;
L = 3;
N = 50;
PdB = [0, -2, -4];
PndB = 10;
[A, S] = generateAVec(M, N, K, PdB, PndB);

R = (1/N)*(A*A');
[U, sval, V] = svd(A);
sval = diag(sval);
[eigvec,eigval0]= eig(R);
[eigval,idx]= sort(diag(eigval0),'descend');
eigvec = eigvec(:,idx);

figure();
stem(sval)
title("Singular Values of A (Part III)")

figure();
stem(eigval)
title("Eigen Values of R (Part III)")

ratio_singVal = sval(3)/sval(4);
ratio_eigVal = eigval(3)/eigval(4);
U_L = U(:, 1:L);
P_S = U_L * U_L.';
P_N = eye(length(P_S), length(P_S)) - P_S;

S_rand = generateRandSVec(100, 20, 20);

S_Music_Known = genSpectrums(S, P_N);
[Music_Known_Stats] = [max(S_Music_Known), mean(S_Music_Known), median(S_Music_Known)];

S_Music_Rand = genSpectrums(S_rand, P_N);
[Music_Rand_Stats] = [max(S_Music_Rand), mean(S_Music_Rand), median(S_Music_Rand)];

% The MUSIC algorithm still works, but the performance is noticeably worse.
% It's now only ~5x greater than the random set.

% Overall, MUSIC is better than MVDR. It contains by far the largest values
% for identifying the source vectors.

%% Part 4 - One Last Thing
newS = S.' * S;

% The Gramian is important because the diagonals tell us the legnth of the
% column vectors, which being 1 signifies that they are normalized (desired
% outcome). The non-diagonal terms show that each component has parts that
% impact the other terms, making it challenging to seperate the sources
% from each other and determine which source constitutes what, unlike if 
% it were to be an orthogonal basis.

function [music] = genSpectrums(S, val)
   for i = 1:size(S,2)
    music(i) = 1/((S(:,i).')*val*S(:,i));
   end
end

function [Sr] = generateRandSVec(M,L,K)
   % Generate S
   Sr = zeros(M, L);
   for i = 1:L
    indicies = randperm(M, K);
    Sr(indicies, i)=abs(randn(1,1));
   end
   Sr = normc(Sr);
end

% PdB is source power
% PndB is noise power
function [A, S] = generateAVec(M, N, K, PdB, PndB)
   dB = 10^(PndB/10);

   % Generate v[n]
   [V] = sqrt(dB) * randn(M, N); % Generate gaussian noise

   % Create L value for "L sources"
   L = length(PdB);
   
   % Generate B
   var_b = 10.^(PdB/10);
   B = zeros(L, N);
   for i = 1:L
    B(i,:) = sqrt(var_b(i))*randn(1, N);
   end

   % Generate S
   S = zeros(M, L);
   for i = 1:L
    indicies = randperm(M, K);
    S(indicies, i) = 1/sqrt(K);
   end

   A = S*B + V.*(1/sqrt(M));
end