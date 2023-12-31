function [optf]=minimax_portfolio()
% 极小极大投资组合优化--
w=0.5; % 权重
xL=0.001; % x的下界

t0=0; % 开始时间
tf=20; % 结束时间
dt=1.5e-3; % 时间步长

% 加载数据并计算协方差矩阵
wk_return=load('wk_return','-ascii'); 
wk_price=load('wk_price','-ascii'); 

week_number = length(wk_return);

% 计算11年的周数
weeks_in_11_years = (11/18)*week_number;
disp([num2str(weeks_in_11_years) ' weeks in 11 years.']);

% year11_weeks = 573; % 2011年的周数

% 截断11年的周回报
trunc_return = wk_return(1:week_number,:);

data_sigma=cov(trunc_return); 
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

% % ---------------------------------------------------------------------

% % 计算投资组合收益的初始化变量
% portfolio_returns = zeros(1, length(tt));
% cumulative_returns = zeros(1, length(tt));

% disp(length(tt)-1)

% % 在循环迭代过程中计算每个周度的投资组合收益，并累积总的周收益：
% for i = 1:length(tt)-1
%     du = NN_MODEL_LIU(xx(:,i), n, xL, P, small_W, W, I_P);
%     xx(:,i+1) = xx(:,i) + (dt) * (du) / 0.001;

%     disp('Size of xx:');
%     disp(size(xx));
%     disp('Size of wk_return:');
%     disp(size(wk_return));

%     % 计算每个周度的投资组合收益
%     portfolio_returns(i) = sum(xx(:, i) .* wk_return(:, i)');
    
%     % 更新投资组合权重
%     nxx(:,i) = FUN_G(xx(:,i), n, xL);
% end

% % 计算最后一个状态的投资组合收益
% portfolio_returns(end) = sum(xx(:, end) .* wk_return(:, end));

% % 计算累积周收益
% cumulative_returns = cumprod(1 + portfolio_returns) - 1;

% count_time=toc; % 计时结束

% % 应用函数FUN_G进行调整
% for i = 1:length(tt)-1
%     xx(:,i)=FUN_G(xx(:,i),n,xL);
% end

% % 对最后一个状态应用FUN_G
% xx(:,length(tt))=FUN_G(xx(:,length(tt)),n,xL);

% % 开始可视化行为
% % ...

%% -----------------------------------------------------------


% 迭代计算
for i = 1:length(tt)-1
    du = NN_MODEL_LIU(xx(:,i),n,xL,P,small_W,W,I_P);
    xx(:,i+1) =xx(:,i)+(dt)*(du)/0.001;
end
count_time=toc

for i = 1:length(tt)-1
    nxx(:,i)=FUN_G(xx(:,i),n,xL);
end


xx(:,length(tt))=FUN_G(xx(:,length(tt)),n,xL);

% 在初始化部分添加以下代码
% 计算初始投资组合价值
initial_weights = xx(n+1:2*n, length(tt)-2); % 使用第一周的权重作为初始权重
initial_prices = wk_price(1, :); % 使用第一周的价格作为初始价格


% 初始化矩阵以存储每周的收益率
weekly_returns = zeros(week_number, 1);

% 初始化投资组合价值
portfolio_value = zeros(week_number, 1);

% 初始化每周的总收益率
total_returns = zeros(week_number, 1);


weights = transpose(xx(n+1:2*n, length(tt) - 1));  % 使用最后计算的权重
next_weights = weights;

% 计算每周的收益率和投资组合价值
for i = 1:week_number


    % error: weights = transpose(xx(n+1:2*n, i));  % 使用第 i 周的权重 (1,58)
    % weights = transpose(xx(n+1:2*n, length(tt) - 1));  % 使用最后计算的权重
    weights = next_weights;


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
        total_returns(i) = (portfolio_value(i) / portfolio_value(1)) - 1;
    end
end



% 计算初始的股票价值
initial_prices = sum(weights .* wk_price(1, :));
% 计算结束的股票价值
end_prices = sum(weights .* wk_price(week_number, :));

% 计算收益率
initial_returns = (end_prices - initial_prices) / initial_prices;

years = 18;
% 根据结束和初始的股票价值，计算平均年化收益率
annualized_returns = (end_prices / initial_prices ) ^ (1/years) - 1;

% 输出
disp('初始的股票价值');
disp(initial_prices);
disp('结束的股票价值');
disp(end_prices);
disp('收益率');
disp(initial_returns);
disp('平均年化收益率');
disp(annualized_returns);


% 可视化每周的收益率
figure;
plot(1:week_number, weekly_returns);
xlabel('周');
ylabel('周收益率');
title('每周收益率');

% 可视化每周的投资组合价值
figure;
plot(1:week_number, portfolio_value);
xlabel('周数');
ylabel('投资组合价值');
title('每周投资组合价值');

% 可视化总收益率
figure;
bar(total_returns);
xlabel('周数');
title('总收益率');

plot_transis = 1;

if plot_transis==1
    % 绘制过渡行为图像
    ttt=(round(count_time/(dt))*(dt))*dt/(length(tt)-1);
    count_tt=0;
    
    for jj=1:length(tt)-1
        count_tt(jj+1)=count_tt(jj)+ttt;
    end
    
    figure
    % 绘制y的过渡行为
    % subplot(1,2,2); % 添加子图
    plot(count_tt(1:length(tt)-1),xx(n+1:2*n,1:length(tt)-1));
    title('Transition Behavior of y');
    ylabel('\it y','FontName','Times New Roman');
    xlabel('time','FontName','Times New Roman');
    xlim([0 0.00020]);
    ylim([0 0.2]);

    % 绘制x的过渡行为
    % subplot(1,2,1); % 添加子图
    % plot(count_tt(1:length(tt)-1),xx(1:n,1:length(tt)-1));
    % title('Transition Behavior of x');
    % ylabel('\it x','FontName','Times New Roman');
    % xlabel('time','FontName','Times New Roman');
    % xlim([0 0.000002]); % 设定为实际的时间范围
    % ylim([xL 1]); % 假设x的上界是1，确保y轴的范围能反映出x的变化

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
