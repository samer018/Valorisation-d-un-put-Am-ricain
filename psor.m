function [x,k,err]=psor(A,B,b,g,x0,eps,kmax)
%- Methode de Newton pour resoudre min(Bx-b,x-g)=0; 

%- Initialisations:
PRINT=0; %- Mettre 1 POUR AFFICHAGES

k=0;
x=x0;
err=eps+1;
n=size(A,1);
hq=b-B*x0;

while( k<kmax & err>eps )
  
  k=k+1;

  xold=x;

x(1)=max(g(1),hq(1)/A(1,1));
  for i=2:1:n
      y=ql(A,hq,x,i);
      x(i)=max(g(i),y);
  end
  hq=b-B*x;
  
 
 
  %- Estimateur pour convergence
  err=norm(min(A*x+B*x-b,x-g),'inf');

  if PRINT
    %- Affichages:
    fprintf('k=%5i, err=%12.6f, n(x-xold)=%10.6f\n',k,err,norm(x-xold)); 
    %scanf('%c');
  end

end