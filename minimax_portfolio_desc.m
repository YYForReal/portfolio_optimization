function [optf]=minimax_portfolio()
% 极小极大投资组合优化--
w=0.5; % 权重
xL=0.01; % x的下界

t0=0; % 开始时间
tf=20; % 结束时间
dt=1.5e-3; % 时间步长

% 加载数据并计算协方差矩阵
wk_return=load('wk_return','-ascii'); 
data_sigma=cov(wk_return); 
n=size(data_sigma,1); % 资产数量
x0=rand(2*(n),1); % 随机初始条件

% 初始化时间和状态变量
tt = t0:dt:tf; 
xx = x0;
B=[zeros(1,n) ones(1,n)]; % 构造矩阵B
P = transpose(B)*inv(B*transpose(B))*B; % 投影矩阵 
small_W=transpose(B)*inv(B*transpose(B)); % 另一种形式的投影矩阵
I = eye(n*2); % 单位矩阵
H=eye(n);
H=H*(w-1); % 耦合矩阵H，实际上是-0.5的单位阵
data_sigma_new=2*w*data_sigma; % 调整后的协方差矩阵
Q=zeros(n); % 零矩阵
W = [Q -H;transpose(H) data_sigma_new]; % 组合矩阵W
I_P=I-P; % 计算I-P
                                                                                                                                                                                                     
% 迭代计算
for i = 1:length(tt)-1
    du = NN_MODEL_LIU(xx(:,i),n,xL,P,small_W,W,I_P);
    xx(:,i+1) =xx(:,i)+(dt)*(du)/0.001;
end
count_time=toc; % 计时结束

% 应用函数FUN_G进行调整
for i = 1:length(tt)-1
    xx(:,i)=FUN_G(xx(:,i),n,xL);
end

% 对最后一个状态应用FUN_G
xx(:,length(tt))=FUN_G(xx(:,length(tt)),n,xL);

% 可视化过渡行为
plot_transis=1;

if plot_transis==1
    % 绘制过渡行为图像
    ttt=(round(count_time/(dt))*(dt))*dt/(length(tt)-1);
    count_tt=0;
    
    for jj=1:length(tt)-1
        count_tt(jj+1)=count_tt(jj)+ttt;
    end
    
    figure
    % 绘制y的过渡行为
    subplot(1,2,2); % 添加子图
    plot(count_tt(1:length(tt)-1),xx(n+1:2*n,1:length(tt)-1));
    title('Transition Behavior of y');
    ylabel('\it y','FontName','Times New Roman');
    xlabel('time','FontName','Times New Roman');
    xlim([0 0.00020]);
    ylim([0 0.2]);

    % 绘制x的过渡行为
    subplot(1,2,1); % 添加子图
    plot(count_tt(1:length(tt)-1),xx(1:n,1:length(tt)-1));
    title('Transition Behavior of x');
    ylabel('\it x','FontName','Times New Roman');
    xlabel('time','FontName','Times New Roman');
    xlim([0 0.000002]); % 设定为实际的时间范围
    ylim([xL 1]); % 假设x的上界是1，确保y轴的范围能反映出x的变化

end


end

% NN_MODEL_LIU函数定义
function du = NN_MODEL_LIU(u,n,xL,P,small_W,W,I_P)
new_u=zeros(n*2,1);
x=u(1:n);
y=u(n+1:2*n);

% 限制x和y的值在规定范围内
x=min(1,max(xL,x));
y=min(1,max(0,y));
new_u(1:n)=x;
new_u(n+1:n*2)=y;

% 计算增量du
du = -P*new_u-(I_P)*(u-new_u+W*((I_P)*new_u+small_W))+small_W;
end

% FUN_G函数定义
function u = FUN_G(u,n,xL)
x=u(1:n);
y=u(n+1:2*n);
xR = 1;
yL = 0;
yR = 1;
% 限制x和y的值在规定范围内
x=min(xR,max(xL,x));
y=min(yR,max(yL,y));
u(n+1:n*2)=y;
u(1:n)=x;
end
