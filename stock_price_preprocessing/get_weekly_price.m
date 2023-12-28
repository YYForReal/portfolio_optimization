function get_weekly_price
% 函数用于计算股票的周度价格和周度收益率

% 加载包含股票数据的文件
% load('stocks_weekly_FTSE.mat');
% load('clean_FTSE100_stocks.mat')
dataset_name = "FTSE100_new";
load('clean_stock_'+ dataset_name + '.mat');

% 初始化变量
wk_price = zeros(size(new_stocks, 1), size(new_stocks, 2));
wk_return = zeros(size(new_stocks, 1), size(new_stocks, 2));
mean_return = zeros(size(new_stocks, 2), 1);
stdDev_return = zeros(size(new_stocks, 2), 1);

% 循环处理每支股票
for i = 1:size(new_stocks, 2)
    % 提取每支股票的调整后收盘价数据
    wk_price(:, i) = new_stocks(i).AdjClose;
    
    % 计算每支股票的周度收益率
    wk_return(:, i) = price2ret(new_stocks(i).AdjClose);
    
    % 计算每支股票的平均收益率和标准差
    mean_return(i, :) = mean(wk_return(:, i));
    stdDev_return(i, :) = std(wk_return(:, i));
end

% 计算股票之间的相关性矩阵
Correlation = corrcoef(wk_return);

% 保存相关数据到文件(ASCII编码)
save(dataset_name + 'correlation_matrix', 'Correlation', '-ASCII');
save(dataset_name + 'mean_return', 'mean_return', '-ASCII');
save(dataset_name + 'stdDev_return', 'stdDev_return', '-ASCII');
save(dataset_name + 'wk_price', 'wk_price', '-ASCII');
save(dataset_name + 'wk_return', 'wk_return', '-ASCII');

% 保存日期信息到ASCII文件
date = new_stocks(1).Date;
dn = datenum(date, 'yyyy-mm-dd');
save('date', 'dn', '-ASCII');

end
