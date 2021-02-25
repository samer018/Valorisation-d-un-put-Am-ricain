function [ind,Pex]=Calculref(w,y,x,v)
%Calculer la valeur exacte de l'option put américaine pour x prix et date t


n=size(x,1);

%j=find(b<c & a>=c);
Pex=0;
d=y(1500);
a=zeros(n,1);
b=zeros(n,1);
%c=zeros(I,1);
for i=1:n
    a(i)=x(i)-v/2;
    b(i)=x(i)+v/2;
end


j=n;
for l=1:n
 if ((a(l)<=d) & (b(l)>=d))
            j=min(j,l);                  %trouver i tq y(320) appartient à [x(i)-v/2,x(i)+v/2]
 end
        
end

%s=find(a<c & b>=c);
Pex=w(1500,1);
ind=j;
end



