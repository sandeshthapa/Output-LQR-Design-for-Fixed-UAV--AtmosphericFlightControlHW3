%% MAE 5010- Atmospheric Flight Control HW3
% LQ Regulator Prob 2
% Sandesh Thapa

clear all;
close all
clc;
%%
A = [0 1;0 0];
B = [0 1]';
C = [1 1];
Q = eye(2);
x0 = [1 0.5]';
X = x0*x0';
R = 5.0;

[K,J,P] = OutputLQRProb2(A,B,C,Q,R,x0);
%%
sim('Sim_OutputLQR')
%%
figure
plot(t,x_lqrA(:,1),'--k',t,x_lqrB(:,1),':b',t,x_lqrC(:,1),'-.r',t,x_lqrD(:,1),':m','Linewidth',1.5);
legend('R = 0.1','R = 0.5','R = 1.0','R = 5.0')
xlabel('time (s)')
ylabel('x1')
hold on 

figure
plot(t,x_lqrA(:,2),'--k',t,x_lqrB(:,2),':b',t,x_lqrC(:,2),'-.r',t,x_lqrD(:,2),':m','Linewidth',1.5);
legend('R = 0.1','R = 0.5','R = 1.0','R = 5.0')
xlabel('time (s)')
ylabel('x2')
hold on 

%% 
figure
plot(t,x_lqrA(:,1));
legend('R = 0.1')
xlabel('time (s)')
ylabel('x1')
hold on 

figure
plot(t,x_lqrA(:,2));
legend('R = 0.1')
xlabel('time (s)')
ylabel('x2')
hold on

