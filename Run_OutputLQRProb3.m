clear all;
close all
clc;
%%
A = [-0.131150  0.14858   0.32434     -0.93964;
          0.0       0.0       1.0      0.33976; 
      -10.614        0.0   -1.1793       1.023; 
      0.99655         0.0   -0.001874  -0.25855];

B = [0.00012    0.00032897;
     0.0        0.0;
     -0.1031578 0.020987; 
     -0.0021330  -0.010715];
   
C = [0 0 57.29578 0;
     0 0 0 57.29578];
 
 D = zeros(2,2);
 
qdr = 50; 
qr = 100;
% v = [qdr,qr,qr,qdr,0,0,0];
Q = diag([qdr,qr,qr,qdr]);
x0 = [1.0 1.0 1.0 1.0]';
X = x0*x0';
rho = 1.0;
n = 3;
% for j = 1:n 
R = rho*eye(2);

[K,J,P] = OutputLQRProb3(A,B,C,Q,R,x0);

%%
% K = K_k;
sim('Sim_OutputLQR')
%%
% end 
% for j = 1:n
figure(1)
plot(t,x_lqrA(:,1),t,x_lqrB(:,1),t,x_lqrC(:,1),'Linewidth',1.0);
legend('\rho = 0.1','\rho = 0.5','\rho = 1.0')
xlabel('time (s)')
ylabel('\beta (rad)')
hold on 

figure(2)
plot(t,x_lqrA(:,2),t,x_lqrB(:,2),t,x_lqrC(:,2),'Linewidth',1.0);
legend('\rho = 0.1','\rho = 0.5','\rho = 1.0')
xlabel('time (s)')
ylabel('\phi (rad)')
hold on 

figure(3)
plot(t,x_lqrA(:,3),t,x_lqrB(:,3),t,x_lqrC(:,3),'Linewidth',1.0);
legend('\rho = 0.1','\rho = 0.5','\rho = 1.0')
xlabel('time (s)')
ylabel('p (rad/s)')
hold on 

figure(4)
plot(t,x_lqrA(:,4),t,x_lqrB(:,4),t,x_lqrC(:,4),'Linewidth',1.0);
legend('\rho = 0.1','\rho = 0.5','\rho = 1.0')
xlabel('time (s)')
ylabel('r (rad/s)')
hold on 

