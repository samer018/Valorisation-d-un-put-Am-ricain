function ploot(t,s,P)
global K r sigma vl vr v0
global Xmin Xmax Ymin Ymax
global Smin Smax



figure(1);
clf;


sgraph  =[Xmin;s;Xmax];
Pgraph  =[vl(t);P; vr(t)];


%- European case:
%plot(sgraph,P0(sgraph),'red--'); % Payoff function
%hold on;
%plot(sgraph,Pexgraph,'black.-'); 
%plot(sgraph,Pgraph,'blue.-');
%- American case:
PAYOFF=plot(sgraph,v0(sgraph),'blue.-','Linewidth',2); % Payoff function
hold on;
NUM=plot(sgraph,Pgraph,'black.-','Linewidth',2); 
legend([PAYOFF,NUM],'Payoff','Scheme','Location','Best');

titre=strcat('t=',num2str(t)); title(titre);
xlabel('x');
ylabel('price');
grid;

