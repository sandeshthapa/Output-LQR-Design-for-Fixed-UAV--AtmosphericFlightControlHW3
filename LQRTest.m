%%
clear all; 
close all; 
clc; 

A = [0 1;0 0];
B = [0 1]';
C = [1 1];
Q = eye(2);
x0 = [0 1]';
X = x0*x0';
R = .11;

n = 10000;

K0 = 0.01;
% Jt = zeros(n,1);
Kt = zeros(n,1);
% del_J = zeros(n,1);
% for m = 1:n

A_c = A - B*K0*C; 

S = Q + C'*K0'*R*K0*C;

P = lyap(A_c,S);

Lamda = lyap(A_c,X);
% X = x0*x0';

J = 0.5*trace(P*X);
% Kt(m) = J(m);

Del_k = inv(R)*B'*P*Lamda*C'*inv(C*Lamda*C') - K0;

for m = 2:1000
%     A_c = A - B*K0*C;
%     
%     S = Q + C'*K0'*R*K0*C;
%     
%     P = lyap(A_c,S);
%     
%     Lamda = lyap(A_c,X);
%     % X = x0*x0';
%     
%     J = 0.5*trace(P*X);
%     % Kt(m) = J(m);
%     
%     Del_k = inv(R)*B'*P*Lamda*C'*inv(C*Lamda*C') - K0;
    alpha = 0.1;
%     K_nn = K0 + alpha*Del_k;
%     Kt(m) = K(m);
    ee = eig(A_c);
%     if m > 1
        J_n = 0.5*trace(P*X);
        del_J = J - J_n;
        if  max(real(ee)) & (del_J) < 0.002
            alpha = 0.01;
            K_new(m) = K0 + alpha*Del_k;
%             Kt(m) = K_nn(m);
            break
        else
            m = m+1;
            A_c = A- B*K0*C;
            S = Q + C'*K0'*R*K0*C;
            
            P = lyap(A_c,S);
            
            X = x0*x0';
            
            Lamda = lyap(A_c,X);
            
            J_n = 0.5*trace(P*X);
            Jt(m) = J(m);
            
            Del_k = inv(R)*B'*P*Lamda*C'*(C*Lamda*C') - K0;
            
            X = x0*x0';
            
            Lamda = lyap(A_c,X);
            
            alpha = alpha/4;
            K_new(m) = K0 + alpha*Del_k
%              Kt(m) = K_nn(m);
        end
%     end
K0 = K_new(m);
J = J_n(m);
end
% J = J_n;

x0 = [0 1]';
K = K_new(end);
sys=ss(A-B*K*C,B,C,0);
[y,t,x]=initial(sys,x0,linspace(0,200,1000));
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
R = 2;

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
% x0 = [1 0.1 0]';
X = x0*x0';

J(m) = 0.5*trace(P*X);
Jt(m) = J(m);

Del_k = inv(R)*B'*P*Lamda*C'*inv(C*Lamda*C') - K0;
alpha = 0.01;

K(m) = K0 + alpha*Del_k;
Kt(m) = K(m);
eA = eig(A_c);
if m > 1
    del_J = (Jt(m) - Jt(m-1));
    if  max(real(eA)) < 0 &&  (del_J) < 0.0000002
         K0 = K(m+1);
         J = J(m+1);
        break
%     end 
     
%     if max(real(eA)) < 0 &&  (del_J) < 0.0000002
%         K0 = K(m);

    else 
        m = m+1;
        A_c = A- B*K0*C; 
        S = Q + C'*K0'*R*K0*C;
       
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
K0 = K(m+1);
J = J(m+1);
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

% figure 
% plot(J)
% xlabel('no of iterations')
% ylabel('J')
% 
% hold on 
% figure
% plot(K)
% xlabel('no of iterations')
% ylabel('Gain K ')


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