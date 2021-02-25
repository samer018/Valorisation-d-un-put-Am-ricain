function [U,L]=uldecomp(A)
%---------------------------------------------------
%- DECOMPOSITION  A = U * L
%- case A= tridiag(.,.,.); 
%- with L=tridiag(l,d,0), U=tridiag(0,1,u)
%---------------------------------------------------
n=size(A,1);
d=zeros(n,1);
l=zeros(n-1,1);
u=zeros(n-1,1);
d(n)=(A(n,n)); 
for i=n-1:-1:1
  l(i)=(A(i+1,i));
  u(i)=(A(i,i+1)/d(i+1));
  d(i)=(A(i,i))-l(i)*u(i); 
end
L=diag(d)+diag(l,-1);
U=eye(size(A))+diag(u,1);

