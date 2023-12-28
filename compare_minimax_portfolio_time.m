function [] = compare_minimax_portfolio()
    % 创建一个新的图像
    figure;
    
    % 加载数据并计算协方差矩阵
    wk_return = load('wk_return','-ascii'); 
    wk_price = load('wk_price','-ascii'); 
    
    % 测试不同数据截断情况的执行时间
    % data_sizes = [200, 400, 600, 800];
    % data_sizes = 100:100:900;
    data_sizes = [573,930];
    num_trials = 100; % 外层循环次数

    execution_times = zeros(num_trials, length(data_sizes));
    final_total_returns = zeros(num_trials, length(data_sizes));

    for trial = 1:num_trials
        disp(['Trial ' num2str(trial)]);
        for i = 1:length(data_sizes)
            data_size = data_sizes(i);
            disp("data_size")
            disp(data_size)

            % 截断数据
            wk_return_truncated = wk_return(1:data_size, :);
            wk_price_truncated = wk_price(1:data_size, :);

            % 测试执行时间
            tic;
            [~, ~, total_returns] = minimax_portfolio(wk_return_truncated, wk_price,wk_return);
            execution_times(trial, i) = toc;
            disp(['Execution time for data size ' num2str(data_size) ': ' num2str(execution_times(trial, i)) 's']);

            % 输出最后的total_return
            % 记录最后的total_return
            final_total_returns(trial, i) = total_returns(end);
            disp(['最后的总收益率: ' num2str(final_total_returns(trial, i))])


        end
    end

    % 计算平均执行时间
    average_execution_times = mean(execution_times);

    % 可视化平均执行时间
    plot(data_sizes, average_execution_times, '-o');
    xlabel('数据大小');
    ylabel('平均执行时间 (s)');
    title('平均执行时间对比');

    % 计算平均最后总收益率
    average_final_total_returns = mean(final_total_returns);

    % 可视化平均最后总收益率
    figure;
    plot(data_sizes, average_final_total_returns, '-o');
    xlabel('数据大小');
    ylabel('平均最后总收益率');
    title('平均最后总收益率对比');



end




function [weekly_returns, portfolio_value, total_returns] = minimax_portfolio(wk_return, wk_price,all_wk_return)
% 极小极大投资组合优化

w = 0.5;
xL = 0.005;
t0 = 0; % 开始时间
tf = 20; % 结束时间
dt = 1.5e-3; % 时间步长

% 加载数据并计算协方差矩阵
% wk_return = load('wk_return','-ascii'); 
% wk_price = load('wk_price','-ascii'); 

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

wk_return = all_wk_return;
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

disp("week_number")
disp(week_number)


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
