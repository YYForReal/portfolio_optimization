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

% % -----------------------------------------------------------

%% 这里是tt时间步长版本2
% Initialize matrix to store the value of each asset over time
asset_values = zeros(n, length(tt));

% Initialize matrix to store the portfolio value over time
portfolio_value = zeros(1, length(tt));

% Calculate the initial portfolio value
initial_portfolio_value = sum(x0(1:n));

% Initialize matrix to store the cumulative returns over time
cumulative_returns = zeros(1, length(tt));



wk_return_first_asset = zeros(1, length(tt));

% Iterative calculation
for i = 1:length(tt)-1
    % Calculate the increment in portfolio weights
    du = NN_MODEL_LIU(xx(:,i), n, xL, P, small_W, W, I_P);
    
    % Update portfolio weights using the calculated increment
    xx(:,i+1) = xx(:,i) + (dt)*(du)/0.001;


    % Ensure that indexing is within bounds for wk_return
    index = min(i, size(wk_return, 1));

    % Calculate the value of each asset at time i
    asset_values(:, i) = xx(1:n, i) .* wk_return(index, :)';
       
    

    % Update portfolio value at time i
    portfolio_value(i) = sum(asset_values(:, i));

    % Calculate cumulative return at time i
    cumulative_returns(i) = (portfolio_value(i) / initial_portfolio_value);


    % Normalize weights for the next iteration (optional, depending on your needs)
    % xx(:, i+1) = FUN_G(xx(:, i+1), n, xL);

    % Save wk_return at each time step
    wk_return_first_asset(i) = wk_return(index, 1);

end

% Ensure that indexing is within bounds for wk_return at the last time step
last_index = min(length(tt), size(wk_return, 1));

% Calculate the value of each asset at the last time step
asset_values(:, end) = xx(1:n, end) .* wk_return(last_index, :)';

% Update portfolio value at the last time step
portfolio_value(end) = sum(asset_values(:, end));

count_time=toc; % Stop the timer



% Visualize the transient behavior if required
plot_transis=1;

if plot_transis==1
    % Plot transient behavior
    ttt=(round(count_time/(dt))*(dt))*dt/(length(tt)-1);
    count_tt=0;

    for jj=1:length(tt)-1
        count_tt(jj+1)=count_tt(jj)+ttt;
    end
    % 绘制资产配比
    figure;
    plot(count_tt(1:length(tt)-1),xx(n+1:2*n,1:length(tt)-1));
    hold on;
    ylabel('\it y','FontName','Times New Roman');
    xlabel('time','FontName','Times New Roman');
    ylim([0 0.2]);


    % Visualize the cumulative returns over time
    figure;
    plot(tt, cumulative_returns);
    xlabel('Time');
    ylabel('Cumulative Returns');
    title('Cumulative Returns Over Time');

    % Visualize how wk_return for the first asset changes over time
    figure;
    plot(tt, wk_return_first_asset);
    xlabel('Time');
    ylabel('wk\_return for the First Asset');
    title('wk\_return for the First Asset Over Time');

    % Visualize the portfolio value over time
    figure;
    plot(tt, portfolio_value);
    xlabel('Time');
    ylabel('Portfolio Value');
    title('Portfolio Value Over Time');
end









% if plot_transis==1
%     % 绘制过渡行为图像
%     ttt=(round(count_time/(dt))*(dt))*dt/(length(tt)-1);
%     count_tt=0;
    
%     for jj=1:length(tt)-1
%         count_tt(jj+1)=count_tt(jj)+ttt;
%     end
    
%     figure
%     % 绘制y的过渡行为
%     % subplot(1,2,2); % 添加子图
%     plot(count_tt(1:length(tt)-1),xx(n+1:2*n,1:length(tt)-1));
%     title('Transition Behavior of y');
%     ylabel('\it y','FontName','Times New Roman');
%     xlabel('time','FontName','Times New Roman');
%     xlim([0 0.00020]);
%     ylim([0 0.2]);

%     % 绘制x的过渡行为
%     % subplot(1,2,1); % 添加子图
%     % plot(count_tt(1:length(tt)-1),xx(1:n,1:length(tt)-1));
%     % title('Transition Behavior of x');
%     % ylabel('\it x','FontName','Times New Roman');
%     % xlabel('time','FontName','Times New Roman');
%     % xlim([0 0.000002]); % 设定为实际的时间范围
%     % ylim([xL 1]); % 假设x的上界是1，确保y轴的范围能反映出x的变化

% end


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
