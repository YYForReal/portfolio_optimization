% load("stock_FTSE100.mat");
% load()

% file_name = "stock_HSI_new.mat";
file_name = "stock_FTSE100_new.mat";
output_file_name = "clean_" + file_name; 
load(file_name)


% 计算最大日期长度
maxDateLength = max(arrayfun(@(x) length(x.Date), stocks));

% 创建一个新结构体数组用于存储筛选后的股票
new_stocks = struct('Ticker', {}, 'Date', {}, 'Open', {}, 'High', {}, 'Low', {}, 'Close', {}, 'AdjClose', {}, 'Volume', {});

% 遍历每个股票，筛选出日期长度等于最大长度的股票
for i = 1:length(stocks)
    ticker = stocks(i).Ticker;
    dates = stocks(i).Date;
    % 检查日期长度和是否包含NaN
    if length(dates) == maxDateLength 
        % 将符合条件的股票添加到新的结构体数组中
        new_stocks(end + 1) = stocks(i);
    end

end

% 保存新的结构体数组到.mat文件
save(output_file_name, 'new_stocks');
