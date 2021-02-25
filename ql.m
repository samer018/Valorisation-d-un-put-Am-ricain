function y=ql(A,b,x,q)
%-- Calculate a sum for PSOR method
y=b(q);
for k=1:q-1
  y=y-A(q,k)*x(k); 
end
y=y/A(q,q);