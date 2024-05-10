clear
clc
%%

n=2^7;
sigma=-1
lambda=2;
alpha=-4;
beta=1;
gamma=1;
A=toeplitz([alpha lambda sigma zeros(1,n-3)],[alpha beta gamma zeros(1,n-3)]);

x_star=ones(n,1);
b=A*x_star;

%%
for i=1:10
    [x_fast,time(i),re_fast]=Pentadiagonal_Toeplitz_Fast_Solver(sigma,lambda,alpha,beta,gamma,b,n);
end
t_ave_fast=sum(time)/10;

%%

for i=1:10
    
    [L,U,P]=lu(A);
    %L=P'*L;
    tic
    y=forward_Substitution_System_Solver(L,b);
    xlu=Backward_Substitution_System_Solver_lu(U,y);
    t(i)=toc;
end
t_ave_lu=sum(t)/10;
re_lu=norm(b-A*xlu)/norm(b);
%%
t_ave_fast
t_ave_lu
re_fast
re_lu


