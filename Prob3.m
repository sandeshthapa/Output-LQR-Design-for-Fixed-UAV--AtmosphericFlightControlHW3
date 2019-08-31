close all;clear all; clc
% Prob 3-HW2- MAE 5010- Atmospheric Flight Controls
% Sandesh Thapa

%%
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
x0 = [1 0.1 0]';
X = x0*x0';

n = 1000;

K0 = [0.1];
Jt = zeros(n,1);
Kt = zeros(n,1);
for m = 1:n

A_c = A- B*K0*C; 

S = Q + C'*K0'*R*K0*C;

P = lyap(A_c,S);

Lamda = lyap(A_c,X);
x0 = [1 0.1 0]';
X = x0*x0';

J(m) = 0.5*trace(P*X);
Jt(m) = J(m);

Del_k = inv(R)*B'*P*Lamda*C'*inv(C*Lamda*C') - K0;
alpha = 0.01;

K(m) = K0 + alpha*Del_k;
Kt(m) = K(m);
eA_c = eig(A_c);
if m > 1
    del_J = (Jt(m) - Jt(m-1)); 
    if abs(del_J) < 0.0000002
        break
    end 
     
    if max(real(eA_c)) < 0 &&  (del_J) < 0.0000002
        K0 = K(m);

    else 
        m = m+1;
        A_c = A- B*K0*C; 
        S = Q + C'*K0'*R*K0*C;
        x0 = [1 0.1 0]';
        
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

sim('Sim_Prob3')

figure
plot(t,x_olqr(:,1),'--k',t,x_olqr(:,2),'b',t,x_olqr(:,3),'-.r');
legend('\Deltap','\Delta\phi','\Delta\delta_a','Location','best','Interpreter','latex')
xlabel('time (s)')
ylabel('States')
hold on 

figure 
plot(t,y_olqr)
xlabel('time (s)')
ylabel('Voltage input to the aileron actuator servomotor')


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

%%