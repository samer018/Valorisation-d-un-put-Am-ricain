function y=descente_p(L,b,g)
%-- solving  Ly-b >=0, y>=g, (Ly-b,y-g)=0 (that is , min(Ly-b, y-g)=0),
%-- in the case L_{ii}>0 and  L up triangular, tridiagonal matrix.
n=length(b);
y=zeros(n,1);
y(1)=b(1)/L(1,1); 
% completer y(1)
y(1)=max(y(1),g(1));
for k=2:n
  y(k)=(b(k)-L(k,k-1)*y(k-1))/L(k,k); 
  % completer y(k)
  y(k)=max(y(k),g(k));
end



