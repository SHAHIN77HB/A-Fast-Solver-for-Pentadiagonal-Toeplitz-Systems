function [x,time,relative_error] = Pentadiagonal_Toeplitz_Fast_Solver (sigma,lambda,alpha,beta,gamma,b,n)

A=toeplitz([alpha lambda sigma zeros(1,n-3)],[alpha beta gamma zeros(1,n-3)]);

J=[zeros(n-2,2) eye(n-2);eye(2) zeros(2,n-2)];
A_bar=J*A;
b1=b(1);
b2=b(2);
b3=b(3:end);
A11=A_bar(1:n-2,1:n-2);
w=A_bar(n-1,1:n-2);
s=A_bar(n,1:n-2);
p=A_bar(1:n-2,n-1);
r=A_bar(1:n-2,n);
%% Solving A11 y = r
tic
u=Backward_Substitution_System_Solver(sigma,lambda,alpha,beta,gamma,b3);
v=Backward_Substitution_System_Solver(sigma,lambda,alpha,beta,gamma,p);
z=Backward_Substitution_System_Solver(sigma,lambda,alpha,beta,gamma,r);

%% Computing x

x3=(s*u-b2-(s*v)/(w*v)*(w*u-b1))/(s*z-(s*v)/(w*v)*w*z);
x2=(w*u-b1)/(w*v)-(w*z)/(w*v)*(s*u-b2-(s*v)/(w*v)*(w*u-b1))/(s*z-(s*v)/(w*v)*w*z);

x1=u-v*x2-z*x3;
time=toc;
x=[x1;x2;x3];
relative_error=norm(b-A*x)/norm(b);