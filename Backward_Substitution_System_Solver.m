function y=Backward_Substitution_System_Solver(sigma,lambda,alpha,beta,gamma,rhs)

sigma1=1/sigma;
alpha1=alpha/sigma;
beta1=beta/sigma;
gamma1=gamma/sigma;
lambda1=lambda/sigma;

n=length(rhs);
y(n)=sigma1*rhs(n);
y(n-1)=sigma1*rhs(n-1)-lambda1*y(n);
y(n-2)=sigma1*rhs(n-2)-lambda1*y(n-1)-alpha1*y(n);
y(n-3)=sigma1*rhs(n-3)-lambda1*y(n-2)-alpha1*y(n-1)-beta1*y(n);

for k=n-4:-1:1
    y(k)=sigma1*rhs(k)-lambda1*y(k+1)-alpha1*y(k+2)-beta1*y(k+3)-gamma1*y(k+4);
end
y=y';