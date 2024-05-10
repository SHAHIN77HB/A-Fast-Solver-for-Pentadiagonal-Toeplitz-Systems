function x=Backward_Substitution_System_Solver(A,b)
% Solving an upper triangular system by back-substitution
% A is an n * n upper triangular matrix
% b is n * 1 right hand side vector

x=zeros(length(b),1);
for j=length(b):-1:1
    if (A(j,j)==0)
        error('Matrix is singular!');
    end
    x(j)=b(j)/A(j,j);
    b(1:j-1)=b(1:j-1)-A(1:j-1,j)*x(j);
end