function out=mo_portfolio()
% 双目标神经动力学投资组合优化
% 使用加权和标量化方法
weight=0.15; % 设定权重
tic; % 开始计时

% 加载数据
wk_return=load('wk_return','-ascii'); % 加载周收益数据
mu1=load('mean_return','-ascii'); % 加载平均收益数据
sigma=cov(wk_return); % 计算收益的协方差矩阵
n=size(sigma,1); % 获取资产数量

% 初始化变量
y0=rand(n,1); % 随机生成初始投资组合权重
lambda0=rand; % 随机生成初始lambda值
l=zeros(n,1); % 初始化l，可能用于约束条件
u=ones(n,1); % 初始化u，可能用于约束条件

f0=[y0;lambda0]; % 合并初始条件

% 设置求解常微分方程的参数
step=0.01; % 设置时间步长
tend=3000; % 设置结束时间
tspan=0:step:tend; % 生成时间跨度
[t,f]=ode23(@(t,f) mv_fun(f,n,sigma,mu1,weight), tspan,f0); % 使用ode23求解常微分方程

% 处理优化结果
[r,c]=size(f); % 获取f的尺寸
aa=f(r,1:n); % 获取最后一行的投资组合权重
aa=max(0,aa); % 权重值不小于0
aa=min(1,aa); % 权重值不大于1
f(r,1:n)=aa; % 更新权重值

% 计算最终输出
out(1,1)=f(r,1:n)*mu1; % 计算投资组合的期望收益
out(1,2)=(f(r,1:n)*sigma)*f(r,1:n)'; % 计算投资组合的方差

tt=toc % 结束计时
tt=tt/1000; % 转换计时单位
% 绘制权重随时间变化的图
plot(0:tt/(tend/step):tt,f(:,1:n))
% xlim([0 0.0018])
xlim([0 0.0009]) % yy改

ylim([0 0.2])
xlabel('time') % x轴标签
ylabel('\it y') % y轴标签

% 保存最终权重
final_weight_mo=f(r,1:n)'; % 转置最终的投资组合权重
save('final_weight_mo', 'final_weight_mo','-ASCII'); % 保存到ASCII文件
end
