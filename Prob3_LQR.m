close all;clear all; clc
% Prob 3-HW2- MAE 5010- Atmospheric Flight Controls
% Sandesh Thapa

Lp = -1; % rad/s
Ldela = 30; % s^-2; 
tau = 0.1; % sec 

Deldela_max = 0.436; % rad 
Delphi_max = 0.787; % rad 
delv_max = 10; % volts

A = [Lp 0 Ldela; 
     1  0 0; 
     0  0 -1/tau];

B = [0 0 1/tau]';
C = [0 1 1];

v = [1;1;1];
Q = diag(v);

R = 1;

syms k 
K = k;
A_c = A-B*K*C 
% 
% S = Q + C'*K'*R*K*C;
% x0 = [1 0.1 0]';
% 
% P = lyap(A_c,S)
% 
% J = x0'*P*x0


eig_A_c = eig(A_c)

sum = eig_A_c(1) + eig_A_c(2) + eig_A_c(3);

fun = @(k)- 10*k - 11 ; 
k0 = 0.1;
k_int = fminsearch(fun,k0)

%%


% syms k0 real 
n = 10;
Km = zeros(10,3);
Km(1,:) = [0.1 1 10];
A_c = cell(n,1);
P = cell(n,1);
Jm = zeros(n,1)
for m = 1:n
C_m = diag(C);
A_c{m} = A- B*Km(m,:)*C_m; 
% A_m(1:m+2, 1:m+2) = A_c;

S = Q + C_m'*Km(m,:)'*R*Km(m,:)*C_m;
x0 = [1 0.1 0]';

P{m} = lyap(A_c{m},S)

J(m) = x0'*P{m}*x0;

X = x0*x0';

Lamda = lyap(A_c{m},X)

Im = 0.5*trace(P{m}*X)

del_k = inv(R)*B'*P{m}*Lamda*C'*(C*Lamda*C') - Km(m)

alpha = 0.2;
% eig(A_c(m)
Km = Km + alpha*del_k
JDel = Jmintest(P, X,m,n,J)
end 

% % fun = @(k)A-B*K*C;
% % 
% k = fminsearch(J,x0)

