function out=mo_portfolio()
% Bi-objective neurodynamic portfolio optimization
% Weighted-sum scalarization is used
weight=0.15;
tic;


wk_return=load('wk_return','-ascii');
mu1=load('mean_return','-ascii');
sigma=cov(wk_return);
n=size(sigma,1);

y0=rand(n,1);
lambda0=rand;
l=zeros(n,1);
u=ones(n,1);

f0=[y0;lambda0];

step=0.01;
tend=3000;
tspan=0:step:tend;
[t,f]=ode23(@(t,f) mv_fun(f,n,sigma,mu1,weight), tspan,f0);



[r,c]=size(f);
aa=f(r,1:n);
aa=max(0,aa);
aa=min(1,aa);
f(r,1:n)=aa;

out(1,1)=f(r,1:n)*mu1;
out(1,2)=(f(r,1:n)*sigma)*f(r,1:n)';

tt=toc
tt=tt/1000;
plot(0:tt/(tend/step):tt,f(:,1:n))
xlim([0 0.0018])
ylim([0 0.2])
xlabel('time')
ylabel('\it y')

final_weight_mo=f(r,1:n)';
save('final_weight_mo', 'final_weight_mo','-ASCII');
end