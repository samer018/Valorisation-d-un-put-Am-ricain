function [x,k,err]=newton(B,b,g,x0,eps,kmax)
%- Methode de Newton pour resoudre min(Bx-b,x-g)=0; 

%- Initialisations:
PRINT=0; %- Mettre 1 POUR AFFICHAGES

k=0;
x=x0;
err=eps+1;


while( k<kmax & err>eps )
  
  k=k+1;

  xold=x;

  %- Definition de F(x) et de F'(x):
  F=min(B*x-b,x-g);
  Fp=eye(size(B)); i=find(B*x-b<=x-g); Fp(i,:)=B(i,:);

  %- Definition nouvel x
  x=x-Fp\F;
 
  %- Estimateur pour convergence
  err=norm(min(B*x-b,x-g),'inf');

  if PRINT
    %- Affichages:
    fprintf('k=%5i, err=%12.6f, n(x-xold)=%10.6f\n',k,err,norm(x-xold)); 
    %scanf('%c');
  end

end

