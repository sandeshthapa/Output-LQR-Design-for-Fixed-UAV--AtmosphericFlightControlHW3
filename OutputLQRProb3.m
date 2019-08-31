function [K_k,J_k,P_k] = OutputLQRProb3(A,B,C,Q,R,x0)
K_k = [-0.5 -1.0;0.5 -1.0];
X = x0*x0';
A_k = A - B*K_k*C;
eig(A_k)

Jt = zeros(1000,1);
JtNew = zeros(1000,1);
Kt = zeros(1000,1);
alpha = 0.5;
for k = 1:1000
    %     K_k = K_0;
    A_k = A - B*K_k*C;
    P_k = lyap(A_k,Q + C'*K_k'*R*K_k*C);
    S_k = lyap(A_k,X);
    J_k = 1/2*trace(P_k*X);
    Jt(k) = J_k;
    DeltaK = inv(R)*B'*P_k*S_k*C'*inv(C*S_k*C') - K_k;
    
    for i = 1:1000
        K_kNew = K_k + alpha*DeltaK;
        A_k = A - B*K_kNew*C;
        P_kNew = lyap(A_k,Q + C'*K_kNew'*R*K_kNew*C);
        S_k = lyap(A_k,X);
        EigAk = eig(A_k);
        
        J_kNew = 1/2*trace(P_kNew*X);
        JtNew(i) = J_kNew;
        %       if  i > 1
        %         Delta_J(i) = Jt(i) - Jt(i+1);
        Delta_J = J_k - J_kNew;
        
        if max(real(EigAk)) < 0 && J_kNew < J_k  
            K_k = K_kNew;
%             Kt(k) = K_k;
            J_k = J_kNew;
            Jt(i) = J_k;
            break
        else
            alpha = alpha/10;
            K_kNew = K_k + alpha*DeltaK;
            K_k = K_kNew;
%             Kt(k) = K_k;
            A_k = A - B*K_k*C;
            EigAk = eig(A_k);
            P_kNew = lyap(A_k,Q + C'*K_kNew'*R*K_kNew*C);
            S_k = lyap(A_k,X);
            J_kNew = 1/2*trace(P_kNew*X);
            J_k = J_kNew;
            Jt(i) = J_k;
            DeltaK = inv(R)*B'*P_k*S_k*C'*inv(C*S_k*C') - K_k;
        end
    end
    
    if Delta_J < 0.0000001
        
        K_k = K_kNew;
%         Kt(k) = K_k;
        J_k = J_kNew;
        break
    end
    
end
