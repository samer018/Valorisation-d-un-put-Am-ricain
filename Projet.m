global  K r sigma T Xmin Xmax
K=100; r=0.1; sigma=0.2; T=1; Xmin=0; Xmax=200;

global vl vr v0
v0= @(x) max(K-x,0);    %- Initial values (payoff function) - typically v0=@(s) max(K-s,0)
vl= @(t) K-Xmin;        %- vl= left  value, at Xmin
vr= @(t) 0;             %- vr= right value, at Xmax


%------------------------
%- NUMERICAL DATA
%------------------------
I=2560; N=2560;

% choose SCHEME value in: {'EE-AMER' , 'EI-AMER-UL' , 'EI-AMER-NEWTON', 'EI-AMER-PSOR','EI-AMER-BDF2, 'CN-AMER-NEWTON', CN-AMER-UL'   }
%SCHEME='EE-AMER';
%SCHEME='EI-AMER-UL'; 
%SCHEME='EI-AMER-NEWTON'; 
%SCHEME='EI-AMER-PSOR';
%SCHEME='EI-AMER-BDF2';
%SCHEME='CN-AMER-PSOR';
SCHEME='CN-AMER-NEWTON'
%SCHEME='CN-AMER-UL';
CENTRAGE='CENTRE'; % 'CENTRE', 'DROIT', or 'GAUCHE' 

err_scale=0; %- Echelle pour le graphe d'erreur.
deltan=N/10; %- Eventuellement, Affichage uniquement tous les deltan pas.

%- Affichage des données:
fprintf('sigma=%5.2f, r=%5.2f, Xmax=%5.2f\n',sigma,r,Xmax);
fprintf('Mesh I= %5i, N=%5i\n',I,N);
fprintf('CENTRAGE : %s\n',CENTRAGE)
fprintf('SCHEME: %s\n',SCHEME)




%- DONNEES NUMERIQUES 
%------------------------
%-dt, h (pas d'espace), x (maillage: vecteur colonne), 
dt=T/N;                 %- time step
h=(Xmax-Xmin)/(I+1); 	%- mesh step
x=Xmin+(1:I)'*h;        %- mesh


%- COEFFICIENT CFL
%COMPLETER
cfl=dt/h^2 * (sigma*Xmax)^2;
fprintf('CFL : %5.3f\n',cfl); 
  
  
  
switch CENTRAGE

case 'CENTRE';

  %COMPLETER 
  %alpha(i)=
  %bet(i)=
  %A(i,i), A(i,i-1), A(i,i+1)
  A=zeros(I,I);
  alpha=sigma^2/2 * x.^2 /h^2;
  bet=r*x/(2*h);
  for i=1:I;   A(i,i) = 2*alpha(i) + r; end;
  for i=2:I;   A(i,i-1) = -alpha(i) + bet(i); end;
  for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;

  %q=inline('[(-alpha(1) + bet(1))* ul(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* ur(t)]','t','alpha','bet','I');
  q= @(t) [(-alpha(1) + bet(1))* vl(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* vr(t)];


case 'DROIT' ;

  A=zeros(I,I);
  alpha=sigma^2/2 * x.^2 /h^2;
  bet=r*x/(h);
  for i=1:I;   A(i,i) = 2*alpha(i) + bet(i) + r; end;
  for i=2:I;   A(i,i-1) = -alpha(i) ; end;
  for i=1:I-1; A(i,i+1) = -alpha(i) - bet(i); end;

  q=@(t) [(-alpha(1))* vl(t);  zeros(I-2,1);  (-alpha(end) - bet(end))* vr(t)];

case 'GAUCHE';

  A=zeros(I,I);
  alpha=sigma^2/2 * x.^2 /h^2;
  bet=r*x/(h);
  for i=1:I;   A(i,i)   =2 *alpha(i) - bet(i) + r; end;
  for i=2:I;   A(i,i-1) =  -alpha(i) + bet(i) ; end;
  for i=1:I-1; A(i,i+1) =  -alpha(i) ; end;
  q=@(t) [(-alpha(1)+bet(1))* vl(t);  zeros(I-2,1);  (-alpha(end))* vr(t)];

otherwise 

  fprintf('CENTRAGE not programmed !'); 
end
  
  
%--------------------
%- Initialiser U et graphique
%--------------------
P=v0(x);
ploot(0,x,P);
fprintf('appuyer sur la touche Return'); input('');

%--------------------
%- BOUCLE PRINCIPALE
%--------------------
%- demarrage compteur temps
tic(); 

Id=eye(size(A));


for n=0:N-1

  t=n*dt;

  %- Schema
  switch SCHEME 
 
  case 'EE-AMER'
    % COMPLETE
    P = max(v0(x), (Id - dt*A)*P - dt*q(t));

  case 'EI-AMER-UL';
    if n==0
      B=Id+dt*A; [U,L]=uldecomp(B);
      fprintf('Verification: norm(B-UL)=%10.5f', norm(B-U*L));
      fprintf('\n');
    end
    %- pb: min(Bx-b,x-g)=0, b=Pold-dt*q(t1), g=P0(s);

    t1=t+dt;
    Pold=P-dt*q(t1);
    c=montee(U,Pold);
    P=descente_p(L,c,v0(x));

    %- Verification:
    err=norm(min(B*P-Pold,P-v0(x)));
   % fprintf('Verification: err=%10.5f',err);
    %fprintf('\n');

  case 'EI-AMER-NEWTON';
    if (n==0); B=Id+dt*A; end;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    Pold=P-dt*q(t1);
    %Pold=P;
    %b=P;
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    [P,k]=newton(B,Pold,g,x0,eps,kmax);
    
    %- Verification
    err=norm(min(B*P-Pold,P-v0(x)));
    %fprintf('Verif: err=%10.5f\n',err);
    
  case 'EI-AMER-PSOR';
    if (n==0); B=Id+dt*A; [U,L]=ulldecomp(B); end;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    Pold=P-dt*q(t1);
    %Pold=P; 
    %b=P;
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    [P,k]=psor(L,U,Pold,g,x0,eps,kmax);
    %- Verification
    err=norm(min(B*P-Pold,P-v0(x)));
    %fprintf('Verif: err=%10.5f\n',err);  

   case 'CN-AMER-PSOR';
     if (n==0); B=Id+dt/2*A; [U,L]=ulldecomp(B); end;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    q0=q(t);
    q1=q(t+dt);
    
    Pold=(Id - dt/2*A) * P - dt*(q0+q1)/2;
    %Pold=P;
    %b=P; 
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    [P,k]=psor(L,U,Pold,g,x0,eps,kmax);
    %- Verification
    err=norm(min(B*P-Pold,P-v0(x)));
    %fprintf('Verif: err=%10.5f\n',err);  
    
    case 'EI-AMER-BDF2';
    if (n==0)
    B=Id+dt*A;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    Pold=P-dt*q(t1);
    %Pold=P;
    %b=P;
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    P=newton(B,Pold,g,x0,eps,kmax); 
    Pact=v0(x);
  
    else
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    B=Id+2/3*dt*A;
    Pold=4/3*P-1/3*Pact-2/3*dt*q(t+dt); 
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    Pact=P;
    P=newton(B,Pold,g,x0,eps,kmax);
    end
    
    case 'CN-AMER-NEWTON';
     if (n==0); B=Id+dt/2*A; end;
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    q0=q(t);
    q1=q(t+dt);
    Pold=(Id - dt/2*A) * P - dt*(q0+q1)/2;
    %Pold=P;
    %b=P; 
    x0=P; g=v0(x); eps=1e-10; kmax=50;
    [P,k]=newton(B,Pold,g,x0,eps,kmax);
    %- Verification
    err=norm(min(B*P-Pold,P-v0(x)));
    %fprintf('Verif: err=%10.5f\n',err);  
    
    case 'CN-AMER-UL';
     if (n==0)
     B=Id+dt*A; [U,L]=uldecomp(B); end
    %- pb: min(Bx-b,x-g)=0, b=P, g=P0(s);
    t1=t+dt;
    q0=q(t);
    q1=q(t1);
    Pold=(Id - dt/2*A) * P - dt*(q0+q1)/2;
    c=montee(U,Pold);
    P=descente_p(L,c,v0(x));
    %- Verification
    
    err=norm(min(B*P-Pold,P-v0(x)));
   % fprintf('Verif: err=%10.5f\n',err);  
    
  
  otherwise
    fprintf('SCHEME not programmed'); abort;

  end
  
 
  if mod(n+1,deltan)==0; %- Affichage tous les deltan pas.

   %- Graphiques:
   t1=(n+1)*dt; 
   ploot(t1,x,P); pause(1e-3);
  end
  
end
   [w,k1,h1,y]=Ref(2560,2560);
   %- Calculs d'erreurs:
   [ind,Pex]=Calculref(w,y,x,h);   %- Appel de la solution de reference
   errLI=norm(P(ind)-Pex,'inf');        %- Calcul erreur Linfty
   fprintf('t=%5.2f; Err.Linf=%8.5f',t1,errLI);  
   fprintf('\n');
 
   %input('');
  


  
  
  
  
  
  
  
  