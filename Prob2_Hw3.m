close all;clear all; clc
% Prob 3-HW2- MAE 5010- Atmospheric Flight Controls
% Sandesh Thapa
A = [0 1;0 0];
B = [0 1]';
C = [1 1];
Q = eye(2);

R = 2;
% [K,J,P] = OutputLQR(A,B,C,Q,R);
%%
x0 = [0 1]';
K = K(end);
sys=ss(A-B*K*C,B,C,0);
[y,t,x]=initial(sys,x0,linspace(0,20,1000));
figure(1)
subplot(211)
plot(t,x(:,1)); ylabel('x1')
subplot(212)
plot(t,x(:,2)); ylabel('x2')

%%
A = [0 1;0 0];
B = [0 1]';
C = [1 1];
Q = eye(2);
x0 = [1 1]';
X = x0*x0';
R = 1;

n = 1000;

K0 = 10;
Jt = zeros(n,1);
Kt = zeros(n,1);
del_J = zeros(n,1);
for m = 1:n

A_c = A - B*K0*C; 

S = Q + C'*K0'*R*K0*C;

P = lyap(A_c,S);

Lamda = lyap(A_c,X);
% X = x0*x0';

J(m) = 0.5*trace(P*X);
Jt(m) = J(m);

Del_k = inv(R)*B'*P*Lamda*C'*inv(C*Lamda*C') - K0;

% eA_c = eig(A_c);
% if m > 1
% del_J(m) = (Jt(m) - Jt(m-1)); 
% end 
% 
% if max(real(eA_c)) < 0 & del_J < 0.000002
%     alpha = 0.1; 
%     K(m) = K0 + alpha*Del_k;
%     K = K0 ;
%     J = J(m);
%     Jt(m) = J(m);
%     break 
% else 
%         m = m+1;
%         A_c = A- B*K0*C; 
%         S = Q + C'*K0'*R*K0*C;
%         P = lyap(A_c,S);
%         Lamda = lyap(A_c,X);
%         J(m) = 0.5*trace(P*X);
%         Jt(m) = J(m);
%         
%         Del_k = inv(R)*B'*P*Lamda*C'*(C*Lamda*C') - K0;
%         Lamda = lyap(A_c,X);
%         
%         alpha = alpha/4;
%         K(m) = K0 + alpha*Del_k;
%         K = K0;
% %         J = J(m)
% end
% end 
    
% % 
% % 
alpha = 0.01;
K(m) = K0 + alpha*Del_k;
Kt(m) = K(m);
eA_c = eig(A_c);
if m > 1
    del_J = (Jt(m) - Jt(m-1)); 
    if  max(real(eA_c)) & (del_J) < 0.0000002
        alpha = 0.01;
        K(m) = K0 + alpha*Del_k;
        K0 = K(m);
        break
%     end 
%      
%     if max(real(eA_c)) < 0 &&  (del_J) < 0.0000002
%         K0 = K(m);

    else 
        m = m+1;
        A_c = A- B*K0*C; 
        S = Q + C'*K0'*R*K0*C;
%         x0 = [1 1]';
%         
        P = lyap(A_c,S);
        
        X = x0*x0';

        Lamda = lyap(A_c,X);
        
        J(m) = 0.5*trace(P*X);
        Jt(m) = J(m);
        
        Del_k = inv(R)*B'*P*Lamda*C'*(C*Lamda*C') - K0;
      
        X = x0*x0';
        
        Lamda = lyap(A_c,X);
        
        alpha = alpha/4;
        K(m) = K0 + alpha*Del_k;
        K0 = K(m);
    end
end 
K0 = K(m);
end 

figure 
plot(J)
xlabel('no of iterations')
ylabel('J')

hold on 
figure
plot(K)
xlabel('no of iterations')
ylabel('Gain K ')

%%
x0 = [0 1]';
K = K(end);
sys=ss(A-B*K*C,B,C,0);
[y,t,x]=initial(sys,x0,linspace(0,20,1000));
figure(1)
subplot(211)
plot(t,x(:,1)); ylabel('x1')
subplot(212)
plot(t,x(:,2)); ylabel('x2')


%%
x0 = [0 1]';
K_LQR = K(end);
sim('Sim_Prob3')

figure
plot(t,x_olqr(:,1),'--k',t,x_olqr(:,2),'b');
legend('x1','x2')
xlabel('time (s)')
ylabel('States')
hold on 

% figure 
% plot(t,y_olqr)
% xlabel('time (s)')
% ylabel('Voltage input to the aileron actuator servomotor')


% % fun = @(k)A-B*K*C;
% % 
% k = fminsearch(J,x0)
%%
figure
plot(t,y_olqr,':r',t,V_olqr,'--k','LineWidth',1.5)
legend('Output-LQR', 'LQR')
xlabel('time (s)')
ylabel('Voltage input to the aileron actuator servomotor (V)')
hold on 


figure
plot(t,x_olqr(:,1),'--k',t,y(:,1),'b','LineWidth',1.5);
legend('\Deltap-Output-LQR','\Deltap-LQR','Location','best','Interpreter','latex')
xlabel('time (s)')
ylabel('Roll rate (rad/s)')
hold on 

figure
plot(t,x_olqr(:,2),'--k',t,y(:,2),'b','LineWidth',1.5);
legend('\Delta\phi-Output-LQR','\Delta\phi-LQR','Location','best','Interpreter','latex')
xlabel('time (s)')
ylabel('Roll angle (\phi) (rad/s)','Interpreter','latex')
hold on 

figure
plot(t,x_olqr(:,3),'--k',t,y(:,3),'b','LineWidth',1.5);
legend('\Delta_a-Output-LQR','\Delta_a-LQR','Location','best','Interpreter','latex')
xlabel('time (s)')
ylabel('Aileron deflection angle (rad)')
hold on 