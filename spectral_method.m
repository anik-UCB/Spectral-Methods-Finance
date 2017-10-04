%{ European Put using the spectral methods. This program uses lagdiff.m,...
....lagroots.m, and poldif.m functions to evaluate the differentiation....
    .... matrix, laguerre roots. The program also compares to the classical...
    .... Black Scholes formula
.... @ Aniruddha Dutta, 2017 
%}

%{ -------- Parameters -------------%}

clear all;
clc;
K=40;     % Strike price
v=0.02;   % volatility
r=0.06;   % rate
T=1;      % time to maturity
nt=50;   % time intervals
N=50;     % number of nodes
M=2;      %  number of derivatives required
b=1;      % scaling parameter in Laguerre differentiation

% The function [x, DM] = lagdif(N, M, b) computes the
%  differentiation matrices D1, D2, ..., DM on Laguerre points.
%
[x, DM] = lagdif(N, M, b); 

D1=DM(:,:,1);  % 1-st derivative matrix
D2=DM(:,:,2);  % 2-nd derivative matrix

i=2:N-1; 
I=speye(N);          % unit matrix
X=spdiags(x,0,N,N);  % matrix with x values
X2=X.^2;             % x^2 matrix
H=r*I-r*X*D1-0.5*v*X2*D2; % H operator from Black-Shcholes equation
A=H(i,i);            
xi=x(i);             % price values
yi=max(K-xi,0);      %  Put for European 
h=-T/nt;             % time step
t=T:h:0;             % backward scheme for time (from T to 0)
b=H(i,N)*K;          % Boundary condition
[l,u]=lu(I(i,i)-0.5*h*A); % LU decomposition

for it=2:length(t) 
    % solving Matrix equation using LU decomposition 
    % on every time step
b1=H(i,1)*K*exp(-r*(T-t(it)));      
yi=u\(l\(yi+0.5*h*(A*yi+b+b1)));
b=b1;
end

figure(1)
plot(xi,yi,'-o')
xlabel('Spot Price','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Option Value','FontSize',12,'FontWeight','bold','Color','k')
title('European Put Option vs Spot Price','Fontsize',12,.....
    'Fontweight','bold','Color','k');
legend('Spectral Method');

%{To compare the spectral methods with analytical formula of Black Scholes..
... we have computed the Black Scholes value functions at each step}%


dt = T/N;           % value of time step
t = 0:dt:T;         % time vector
Put1=[];
for i = 1:1:48
    Put_Col = []; 
    [Call,Put] = blsprice(xi(i),40,.06,T-i*dt,.2);  % BS values
    Put_Col = [Put_Col; Put];
    Put1 = [Put1 Put_Col];
    
end

%{
Black Scholes European Put plot

print Put1   % Black Scholes value functions
plot(xi,Put1,'-*')
xlabel('Spot Price','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Option Value','FontSize',12,'FontWeight','bold','Color','k')
title('European Put Option vs Spot Price','Fontsize',12,.....
    'Fontweight','bold','Color','k');
legend('Black Scholes');

%}


