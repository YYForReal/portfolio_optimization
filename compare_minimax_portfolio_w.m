function [] = compare_minimax_portfolio()
% 对比不同 w 超参数的实验
ws = [0.2, 0.5, 0.8, 1.0]; % 不同的 w 超参数

% 创建一个新的图像
figure;

% 执行极小极大投资组合优化和可视化
for i = 1:length(ws)
    w = ws(i);

    % 执行极小极大投资组合优化
    [weekly_returns, portfolio_value, total_returns] = minimax_portfolio(w);

    % 在同一子图中叠加绘图
    subplot(2, 2, 1); % 第一行第一列
    plot(1:length(weekly_returns), weekly_returns, 'o-', 'DisplayName', ['w = ' num2str(w)]);
    hold on;

    subplot(2, 2, 2); % 第一行第二列
    plot(1:length(portfolio_value), portfolio_value, 'o-', 'DisplayName', ['w = ' num2str(w)]);
    hold on;

    % 在新的子图中绘制总收益率
    subplot(2, 2, 3); % 第二行第一列
    plot(1:length(total_returns), total_returns, 'o-', 'DisplayName', ['w = ' num2str(w)]);
    hold on;




end

% 添加图例
subplot(2, 2, 1);
xlabel('周');
ylabel('周收益率');
title('每周收益率');
legend('Location', 'Best');

subplot(2, 2, 2);
xlabel('周');
ylabel('投资组合价值');
title('每周投资组合价值');
legend('Location', 'Best');

% 在总收益率图中添加标题和标签
subplot(2, 2, 3);
xlabel('周');
ylabel('总收益率');
title('总收益率对比');
legend('Location', 'Best');


end

function [weekly_returns, portfolio_value, total_returns] = minimax_portfolio(w)
% 极小极大投资组合优化

xL = 0.01; % x的下界
t0 = 0; % 开始时间
tf = 20; % 结束时间
dt = 1.5e-3; % 时间步长

% 加载数据并计算协方差矩阵
wk_return = load('wk_return','-ascii'); 
wk_price = load('wk_price','-ascii'); 

data_sigma = cov(wk_return); 
n = size(data_sigma,1); % 资产数量
x0 = rand(2*(n),1); % 随机初始条件

% 初始化时间和状态变量
tt = t0:dt:tf;

xx = x0;
B = [zeros(1,n) ones(1,n)]; % 构造矩阵B
P = transpose(B)*inv(B*transpose(B))*B; % 投影矩阵 
small_W = transpose(B)*inv(B*transpose(B)); % 另一种形式的投影矩阵
I = eye(n*2); % 单位矩阵
H = eye(n);
H = H*(w-1); % 耦合矩阵H，实际上是-0.5的单位阵
data_sigma_new = 2*w*data_sigma; % 调整后的协方差矩阵
Q = zeros(n); % 零矩阵
W = [Q -H;transpose(H) data_sigma_new]; % 组合矩阵W
I_P = I - P; % 计算I-P

% 迭代计算
for i = 1:length(tt)-1
    du = NN_MODEL_LIU(xx(:,i), n, xL, P, small_W, W, I_P);
    xx(:,i+1) = xx(:,i) + (dt)*(du)/0.001;
    
    % 对于下一次迭代，可以选择规范化权重
    % nxx(:, i+1) = FUN_G(xx(:, i+1), n, xL);
end

week_number = length(wk_return);

% 计算初始投资组合价值
initial_weights = xx(n+1:2*n, length(tt)-2); % 使用第一周的权重作为初始权重
initial_prices = wk_price(1, :); % 使用第一周的价格作为初始价格

% 初始化矩阵以存储每周的收益率
weekly_returns = zeros(week_number, 1);

% 初始化投资组合价值
portfolio_value = zeros(week_number, 1);

% 初始化每周的总收益率
total_returns = zeros(week_number, 1);

% 计算每周的收益率和投资组合价值
for i = 1:week_number
    % 获取第 i 周的权重
    weights = transpose(xx(n+1:2*n, i));  % 使用第 i 周的权重 (1,58)

    % 计算每支股票在第 i 周的收益率
    stock_returns = (wk_return(i, :)) .* (weights)  ; % (1,58)

    % 计算整个投资组合在第 i 周的收益率
    weekly_returns(i) = sum(stock_returns);

    % 计算第 i 周的投资组合价值
    portfolio_value(i) = sum(weights .* wk_price(i, :));

    % 计算累积总收益率
    if i == 1
        total_returns(i) = 0;
    else
        % total_returns(i) = (1 + total_returns(i-1)) * (1 + weekly_returns(i)) - 1;
        total_returns(i) = portfolio_value(i) / portfolio_value(1);
    end
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
