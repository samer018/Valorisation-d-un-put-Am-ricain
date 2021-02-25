function y=P0(s)
%- payoff function
global K
payoff=1;
switch payoff
case 1; 
  y=max(K-s,0);
case 2; 
  % TO COMPLETE 
  y=zeros(size(s)); i=find(s>=K/2 & s<=K); y(i)=K/2; 
end


