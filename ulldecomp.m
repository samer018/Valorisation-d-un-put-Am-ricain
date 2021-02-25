function [U,L]=ulldecomp(A)
%---------------------------------------------------
%- DECOMPOSITION  A = U * L
%- case A= tridiag(.,.,.); 
%- with L=tridiag(l,d,0), U=tridiag(0,1,u)
%---------------------------------------------------
n=size(A,1);
L=zeros(n,n);
U=zeros(n,n); 
U(1,2)=A(1,2);
L(n,n-1)=A(n,n-1);
L(n,n)=A(n,n);
L(1,1)=A(1,1);
for i=2:1:n-1
U(i,i+1)=A(i,i+1);
L(i,i)=A(i,i);
L(i,i-1)=A(i,i-1);
end
