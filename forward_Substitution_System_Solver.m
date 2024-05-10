function[x]=forward_Substitution_System_Solver(L,b)
%This function solve a lower triangular system using forkward substitution
%method. The standard call is: "x=backward(U,b) in wich L and b represent
%respectively the lower triangular matrix and  the known term.

[m,n]=size(L);
if n~=m
    error('matrix mast be square')
end
x(1)=b(1)/L(1,1);
%bacward substitution
for k=2:m
        x(k)=(b(k)-L(k,1:k-1)*x(1:k-1)')/L(k,k);
end
x=x';
end