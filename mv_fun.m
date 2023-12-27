function df=mv_fun(f,n,sigma,mu1,weight)
beta=10;
e1=10;
e2=10;

y=f(1:n);
lamada=f(n+1);

dfy=2*(1-weight)*(sigma*y)-weight*mu1;
df=dfy;

%equality constraints
h1=sum(y)-1;
h=h1;

%dh1=[ones(n,1);zeros(n,1)];
dh1=ones(n,1);
dh=dh1;

gammamu=diag(lamada.^2);
p=y-(df+dh*lamada+beta*dh*gammamu*h);
%py=min(u,max(l,p(1:n)));
 py=p(1:n);
 py=max(0,py);
 py=min(1,py);

dy=-y+py;
ph=lamada+h;
ph=min(1,ph);
ph=max(0,ph);
df=[e1*dy;e2*(-lamada+ph)];
%df=[e1*dy;e2*(-h)];
end