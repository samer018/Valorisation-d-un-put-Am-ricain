function x=montee(U,y)
%-- solving Ux=y: "montee", with U(k,k)=1 for all k
n=length(y);
x=zeros(n,1); 
x(n)=y(n);
for k=n-1:-1:1
  x(k)=(y(k)-U(k,k+1)*x(k+1)); 
end
