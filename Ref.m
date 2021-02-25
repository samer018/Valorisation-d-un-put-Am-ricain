%-----------------------------------------------
%- Solution de reference 
%-----------------------------------------------

function [w,k,h,y]=Ref(I,N)
global r sigma T Xmax Xmin
global vl vr v0

k=T/N;                 %- time step
h=(Xmax-Xmin)/(I+1); 	%- mesh step
y=Xmin+(1:I)'*h;        %- mesh

A=zeros(I,I);
  alpha=sigma^2/2 * y.^2 /h^2;
  bet=r*y/(2*h);
  for i=1:I;   A(i,i) = 2*alpha(i) + r; end;
  for i=2:I;   A(i,i-1) = -alpha(i) + bet(i); end;
  for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;

  %q=inline('[(-alpha(1) + bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)]','t','alpha','bet','I');
  q= @(t) [(-alpha(1) + bet(1))* vl(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* vr(t)];
  
%- Initialiser U et graphique
%--------------------
P=v0(y);

%--------------------
%- BOUCLE PRINCIPALE

tic(); 

Id=eye(size(A));
w=zeros(I,1);
%Applcation de shéma BDF2 pour approcher la solution exacte.

  for n=0:N-1

  t=n*k;
  if (n==0)
    B=Id+k*A;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+k;
    Pold=P-k*q(t1);
    %Pold=P;
    %b=P;
    x0=P; g=v0(y); eps=1e-10; kmax=50;
    P=newton(B,Pold,g,x0,eps,kmax); 
    Pact=v0(y);
  
  else
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    B=Id+2/3*k*A;
    Pold=4/3*P-1/3*Pact-2/3*k*q(t+k); 
    x0=P; g=v0(y); eps=1e-10; kmax=50;
    Pact=P;
    P=newton(B,Pold,g,x0,eps,kmax);
        
  end
  
  end
  w=P;
