function [optf]=minimax_portfolio()
%minimax portfolio optimzation
w=0.5;
xL=0.01;

t0=0;
tf=20;
dt=1.5e-3;

wk_return=load('wk_return','-ascii');
data_sigma=cov(wk_return);
n=size(data_sigma,1);
x0=rand(2*(n),1);


tt = t0:dt:tf;
xx = x0;
B=[zeros(1,n) ones(1,n)];
P = transpose(B)*inv(B*transpose(B))*B;
small_W=transpose(B)*inv(B*transpose(B));
I = eye(n*2);
H=eye(n);
H=H*(w-1);
data_sigma_new=2*w*data_sigma;
Q=zeros(n);
W = [Q -H;transpose(H) data_sigma_new];
I_P=I-P;

tic


for i = 1:length(tt)-1
    du = NN_MODEL_LIU(xx(:,i),n,xL,P,small_W,W,I_P);
    xx(:,i+1) =xx(:,i)+(dt)*(du)/0.001;
end
count_time=toc

for i = 1:length(tt)-1
    nxx(:,i)=FUN_G(xx(:,i),n,xL);
end


xx(:,length(tt))=FUN_G(xx(:,length(tt)),n,xL);


plot_transis=1;

if plot_transis==1
    
    % plot transient behavior
    ttt=(round(count_time/(dt))*(dt))*dt/(length(tt)-1);
    count_tt=0;
    
    for jj=1:length(tt)-1
        count_tt(jj+1)=count_tt(jj)+ttt;
    end
    
    figure
    plot(count_tt(1:length(tt)-1),xx(n+1:2*n,1:length(tt)-1));
    hold on;
    ylabel('\it y','FontName','Times New Roman');
    xlabel('time','FontName','Times New Roman');
    ylim([0 0.2])
end


end


function du = NN_MODEL_LIU(u,n,xL,P,small_W,W,I_P)
new_u=zeros(n*2,1);
x=u(1:n);
y=u(n+1:2*n);

x=min(1,max(xL,x));
y=min(1,max(0,y));
new_u(1:n)=x;
new_u(n+1:n*2)=y;

du = -P*new_u-(I_P)*(u-new_u+W*((I_P)*new_u+small_W))+small_W;
end

function u = FUN_G(u,n,xL)
x=u(1:n);
y=u(n+1:2*n);
xR = 1;
yL = 0;
yR = 1;
x=min(xR,max(xL,x));
y=min(yR,max(yL,y));
u(n+1:n*2)=y;
u(1:n)=x;
end

