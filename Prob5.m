%% Problem 5(a) 

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
K1 = [  -1.1914   -0.7028;
         0.6179   -1.0284];
     
K2 = [ -1.4810   -0.6327;
    0.6555   -1.0136];

K3 =   [-1.4414   -0.6528;
    0.6468   -1.0087];

G = tf(ss(A,B,C,D));
% G0 = ss(A,B,C,D,2);
% Lp = loopsens(G0,K1)
% figure
% sigma(Lp.To)
% hold on 

figure
sigma(G*K1)
hold on 

figure 
sigma(G*K2);
hold on 

figure 
sigma(G*K3);
hold on 

%% 5(b)
Lp = loopsens(G,K2);

Wii = makeweight(0.30,30,10);
WII = [Wii,Wii;
       Wii,Wii];
   
M = -WII*Lp.To; 
figure
sigma(M)
hold on 

HinfNorm2 = hinfnorm(M)

if HinfNorm2 > 1.0 
    fprintf('Controller is not robust')
else 
    fprintf('System is robutst')
end

%%
Wii = makeweight(0.56,30,10);
WII = [Wii,Wii;
       Wii,Wii];
   
M = -WII*Lp.To; 

HinfNorm2 = hinfnorm(M)

if HinfNorm2 > 1.0 
    fprintf('Controller is not robust')
else 
    fprintf('System is robutst')
end