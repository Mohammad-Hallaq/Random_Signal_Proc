% SECTION A 
% FIRST CASE 
%3
clc
close all
clear

A = 2 ;
sig2 = 0.5;
M = 10;
w = sqrt(sig2)*randn(1,M);
x = A + w;
A_MLE = mean(x);
display(A_MLE);
%%

%4
clc
close all
clear

A = 2 ;
sig2 = 0.5;
M = 100;
k = 500;
x = zeros(k , M);
A_MLE = zeros(1,k);
for i = 1 : k
    w = sqrt(sig2)*randn(1,M);
    x(i,:) = A + w;
    A_MLE(i) = mean(x(i,:));
end
hist(A_MLE);
%%

%5
clc
close all 
clear

M = 3:3:200;
A =2;
sig2 = 0.5;
k = 50;
l = 1;
for i = 3:3:200
    
    for j = 1 :k
   
    w = sqrt(sig2)*randn(1,i);
    x = zeros(k , i);
    
    x(j,:) = A + w;
    
    A_MLE(j) = mean(x(j,1:i));
    
    end
    u(l) = mean(A_MLE);
    S(l) = var(A_MLE);
    l = l+1;
end
plot(M,u );
hold on ;
plot(M ,A*ones(1,length(M)));
figure();
semilogy(M,S, M , sig2./M);



%%
% SECTION A
% SECOND CASE 
%3
clc
close all
clear

A = 2 ;
sig2 = 0.5;
M = 10;
w = sqrt(sig2)*randn(1,M);
x = A + w;
sig2_MLE = (sum((x - A).^2))/M;
display(sig2_MLE);
%%

%4
clc
close all
clear

A = 2 ;
sig2 = 0.5;
M = 100;
k = 500;
x = zeros(k , M);
sig2_MLE = zeros(1,k);
for i = 1 : k
    w = sqrt(sig2)*randn(1,M);
    x(i,:) = A + w;
    sig2_MLE(i) = (sum((x(i,:) - A).^2))/M;
end
hist(sig2_MLE);
%%

%5
clc
close all
clear 

M = 3:3:200;
A = 2;
sig2 = 0.5;
k = 50;
l = 1;
for i = 3:3:200
    
    for j = 1 :k
   
    w = sqrt(sig2)*randn(1,i);
    x = zeros(k , i);
    
    x(j,:) = A + w;
    
    sig2_MLE(j) =  (sum((x(j,1:i) - A).^2))./i;
    
    end
    u(l) = mean(sig2_MLE);
    S(l) = var(sig2_MLE);
    l = l+1;
end
plot(M,u );
hold on ;
plot(M ,sig2*ones(1,length(M)));
figure();
semilogy(M,S, M , sig2./M);

%%
% SECTION A
% THIRD CASE 
%8
clc
close all
clear

A = 2;
sig2 = 0.5;
M = 3:3:2000;
k =1;
for i = 3:3:2000
    w = sqrt(sig2)*randn(1,i);
    x = A + w;
    A_MLE(k) = mean(x);
    sig2_MLE(k) = (sum((x - A_MLE(k)).^2))/i;
    k = k+1;
end
ERR_A = A_MLE - A;
ERR_sig2 = sig2_MLE - sig2;
plot(M,ERR_sig2);
%%

%9
clc
close all
clear

A = 0:0.1:3;
sig2 = 0:0.1:3;
M = length(A);
x = A +randn(1,M).*sqrt(sig2);
p = ones(length(sig2),length(A));
for i=1:length(sig2)
    for j =1:length(A)
    p(i,j) = (1./((2*pi*sig2(i)).^(M/2)).*exp((-1./(2*sig2(i)))*sum((x - A(j)).^2)));
    end
    
end
A_MLE = mean(x);
sig2_MLE = (sum((x - A_MLE).^2))/M;
mesh(sig2,A,p);
title('likelihood function');
ylabel('A');
xlabel('variance');
%%

%SECTION B N(A,A) .................
%3
clear
clc

A = 2;
M = 10;
w = sqrt(A)*randn(1,M);
x = A + w;
A_ML = (-1/2)+sqrt(((sum(x.^2))/M)+ 1/4);
disp(A_ML);

%4

k = 1e4;
A_ML_2 = zeros(1,k);
for i = 1 : 1e4
    w = sqrt(A)*randn(1,M);
    x = A + w;
    A_ML_2(i) = (-1/2)+sqrt(((sum(x.^2))/M)+ 1/4);
end
mean(A_ML_2)

%5

A_ML_3 = zeros(1,k);
j = 1;
for N = 10:10:200
    
    for i = 1 : 1e4
    w = sqrt(A)*randn(1,N);
    x = A + w;
    A_ML_3(i) = (-1/2)+sqrt(((sum(x.^2))/N)+ 1/4);
    end
    Mean_Esm(j) = mean(A_ML_3);
    
    j = j+1;
end

%6

Est_Error = (A_ML_3 - A);
Bias_Error = Mean_Esm - A;
plot(Est_Error);
title("ESTIMATION ERROR");
figure;
plot(10:10:200,Bias_Error);
title("BIAS ERROR");
%%

%SECTION C .................
%2
clc
clear 
close all

f0 = 0.14;
A = 1;
N = 100;
phi = 0.4;
sig2 = 5;
n = 0:N-1;
w = sqrt(sig2)*randn(1,N);
x = A*cos(2*pi*f0*n + phi) + w;
phi_ML = -1*atan(sum(x.*sin(2*pi*f0*n))/sum(x.*cos(2*pi*f0*n)));
x_est =  A*cos(2*pi*f0*n + phi_ML) + w;
figure;
plot(n,x);
hold on
plot(n,x_est);
legend('THE ORIGINAL SIGNAL','THE ESTIMATED SIGNAL');
%%

%3 and 4
clc

N=10:10:200;
SQUARE_ERR = zeros(1,length(N));
mean_phi = zeros(1,length(N));
k =1;
for j=10:10:200
    w = sqrt(sig2)*randn(1,j);
    n = 0:j-1;
    x = A*cos(2*pi*f0*n + phi) + w;
    phi_ML = -atan(sum(x.*(sin(2*pi*f0*n)))/sum(x.*(cos(2*pi*f0*n))));
    SQUARE_ERR(k) = (phi - phi_ML)^2;
    mean_phi(k) = mean(phi_ML);
    k = k+1;
end
BIAS_ERR = mean_phi - phi;
figure;
plot(SQUARE_ERR);
figure;
plot(BIAS_ERR);
%%

%5 and 6 
k = 100;
f0 = 0.14;
A = 1;
N = 30;
phi = 0.4;
sig2 = 5;
n = 0:N-1;
for i = 1 : k
    w = sqrt(sig2)*randn(1,N);
    x = A*cos(2*pi*f0*n + phi) + w;
    phi_ML(i) = -1*atan(sum(x.*sin(2*pi*f0*n))/sum(x.*cos(2*pi*f0*n)));
end
mean_phi = mean(phi_ML)
hist(phi_ML);


