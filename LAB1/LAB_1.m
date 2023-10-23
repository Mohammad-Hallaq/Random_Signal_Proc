clear all;
clc;
close all;

A=1;
sigma=sqrt(0.05);
N = 1000 ;
n = 0 : N-1 ;
CRLB_A = (2*sigma * sigma) ./ n ;
plot( n , CRLB_A);
grid on;
title('CRLB A');
xlabel(' n ');
ylabel('variance');
axis ([0  500  0  0.05] ) ;
%%

clear all;
clc;
close all;

A=1;
sigma=sqrt(0.05);
f0 = 0.08 ;
phi = pi /4 ;
N = 1000 ;
n = 0 : N-1 ;
x=A*cos(2*pi*f0*n+phi)+sqrt(sigma)*randn(1,N);
A_shabo = zeros(1 , 100); 
for J = 1 : N
    A_shabo(1,J) = (2 * sum(x .* cos(2.*pi .* f0 .* n + phi )))./N ;
end
mean_A_shaboo = mean( A_shabo) ;
var_A_shaboo = var(A_shabo);
VARA1 =var(A_shabo - mean(A_shabo)) ;
m1 = A_shabo - A ;
mes_A_shaboo = mean(m1 .* m1);
%%

clear all;
clc;
close all;

A=1;
sigma=sqrt(0.05);
f0 = 0.08 ;
phi = pi /4 ;
i=1;
N=20:20:80;
for J=20:20:80
    n=0:J-1;
    for u=1:1000
    x=A*cos(2*pi*f0*n+phi)+sqrt(sigma)*rand(1,J);
    phi_shabo(u)=-atan(sum(x.*sin(2*pi*f0*n))/sum(x.*cos(2*pi*f0*n)));
    end
    mean_phi_shabo(i)=mean(phi_shabo);
    Var_phi_shabo(i)=var(phi_shabo);
    i=i+1;
end
figure();
plot(N,mean_phi_shabo); 
grid on; 
title('Mean of phi shabo');
figure();
plot(N,Var_phi_shabo); 
grid on; 
title('variance of phi shabo');
%%

clear all;
clc;
close all;

A=1;
sigma=sqrt(0.05);
f0 = 0.08 ;
phi = pi /4 ;
i=1;
N=20:20:5000;
CRLB_phi = (2 * sigma*sigma)./(N * A * A) ;
for J=20:20:5000
    n=0:J-1;
    for u=1:1000
    x=A*cos(2*pi*f0*n+phi)+sqrt(sigma)*rand(1,J);
    phi_shabo(u)=-atan(sum(x.*sin(2*pi*f0*n))/sum(x.*cos(2*pi*f0*n)));
    end
    mean_phi_shabo(i)=mean(phi_shabo);
    Var_phi_shabo(i)=var(phi_shabo);
    i=i+1;
end
figure();
plot(N , phi*ones(1,length(N)) , 'r'); 
hold on ;
plot(N,mean_phi_shabo,'g'); 
grid on; 
title('Mean of phi shabo & the real value of phi');
figure();
plot(N,CRLB_phi , 'r'); 
hold on ;
plot(N,Var_phi_shabo,'g'); 
grid on; 
title('variance of phi shabo & CRLB of phi');
%%

clc;
clear;
close all;

N=80;
f0=0.08;
phi=pi/4;
sigma=sqrt(0.05);
SNR = [ 0 5 10 15] ;
i=1;
A = 10.^(SNR/20) ;
for J=1:length(SNR)
    n=0:N-1;
    for u=1:1000
        x=A(J)*cos(2*pi*f0.*n+phi)+sigma*randn(1,N);
        phi_Shabo(u)=-atan(sum(x.*(sin(2*pi*f0*n)))/sum(x.*(cos(2*pi*f0*n))));
    end
    mean_phi_shabo(i)=mean(phi_Shabo);
    var_phi_shabo(i)=var(phi_Shabo);
    i=i+1;
end
figure();
plot(SNR,mean_phi_shabo);
xlabel('SNR(dB)');
ylabel('practical E(phi)');  
title ('Mean of phi shabo');
grid on;  
figure();
plot(SNR,var_phi_shabo);
xlabel('SNR(dB)');
ylabel('practical var(phi)');
title ('Varince of phi shabo');
grid on;
%%

clc;
clear;
close all;

N=80;
f0=0.08;
phi=pi/4;
sigma=sqrt(0.05);
SNR = [ 0 5 10 15] ;
i=1;
A = 10.^(SNR/20) ;
CRLB_phi = (2 * sigma*sigma)./(N .* A .* A) ;
for J=1:length(SNR)
    n=0:N-1;
    for u=1:1000
        x=A(J)*cos(2*pi*f0.*n+phi)+sigma*randn(1,N);
        phi_Shabo(u)=-atan(sum(x.*(sin(2*pi*f0*n)))/sum(x.*(cos(2*pi*f0*n))));
    end
    mean_phi_shabo(i)=mean(phi_Shabo);
    var_phi_shabo(i)=var(phi_Shabo);
    i=i+1;
end
figure();
semilogy(SNR,mean_phi_shabo,SNR,phi*ones(1,length(SNR)));
xlabel('SNR(dB)');
ylabel('practical E(phi),real value of phi');  
title ('Mean of phi shabo & the real value of phi');
grid on;  
figure();
semilogy(SNR,var_phi_shabo,SNR,CRLB_phi);
xlabel('SNR(dB)');
ylabel('practical var(phi),CRLB phi');
title ('Varince of phi shabo and CRLB of phi');
grid on;
%%

clc;
clear;
close all;

N=80;
n=0:N-1;
f0=0.08;
phi=pi/4;
sigma=sqrt(0.05);
for J=1:1:4
if(J<3) SNR=1;
else SNR=10^(15/20);
end
x= sqrt(2*SNR*sigma*sigma)*cos(2*pi*f0.*n+phi)+sigma*randn(1,N);
phi_shabo=-atan(x.*sum(sin(2*pi*f0*n))/sum(x .*cos(2*pi*f0*n)));
phi2=-pi:0.001:pi;
W=-N/2*log(2*pi*sigma^2)*ones(length(phi2),1)-(1/(2*sigma*sigma)*sum(((x'*ones(1,length(phi2)))-(sqrt(2*SNR*sigma*sigma))*cos(2*pi*f0*(n'*ones(1,length(phi2)))+ones(N,1)*phi2)).^2))';
subplot(2,2,J);
plot(phi2,W);
grid on; 
if(J<3) title('L for SNR=0 dB');
else title('L for SNR=15 dB');
end
xlabel('phi(rad)');
end
