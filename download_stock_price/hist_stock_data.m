function stocks = hist_stock_data(start_date, end_date, varargin)
% 函数用于获取历史股票数据

stocks = struct([]);        % 初始化数据结构

%% 解析输入

% 格式化开始和结束日期为 Posix 时间
origDate = datenum('01-Jan-1970 00:00:00', 'dd-mmm-yyyy HH:MM:SS');
if ischar(start_date)
    startDate = (datenum(start_date, 'ddmmyyyy') - origDate) * 24 * 60 * 60;
else
    startDate = (floor(start_date) - origDate) * 24 * 60 * 60;
end
if ischar(end_date)
    endDate = (datenum(end_date, 'ddmmyyyy') - origDate) * 24 * 60 * 60;
else
    endDate = (floor(end_date) - origDate) * 24 * 60 * 60;
end

% 确定用户是否指定了频率
temp = find(strcmp(varargin,'frequency') == 1);
 if isempty(temp)
    freq = 'd';  % 如果未给出，默认为每日数据
else
    freq = varargin{temp+1};
    varargin(temp:temp+1) = [];
end
clear temp

% 确定用户是否指定了事件类型
temp = find(strcmp(varargin,'type') == 1);
if isempty(temp)
    event = 'history';  % 默认为历史价格
else
    event = varargin{temp+1};
    varargin(temp:temp+1) = [];
end
clear temp

% 如果 varargin 的第一个元素本身是一个单元数组，则假定它是一个股票符号的单元数组
if iscell(varargin{1})
    tickers = varargin{1};
% 否则，检查它是否是 .txt 文件
elseif ~isempty(strfind(varargin{1},'.txt'))
    fid = fopen(varargin{1}, 'r');
    tickers = textscan(fid, '%s'); tickers = tickers{:};
    fclose(fid);
% 否则，假定它要么是单个股票符号，要么是股票符号列表
else
    tickers = varargin;
end

%% 获取历史数据

h = waitbar(0, '请等待...');  % 创建等待条
idx = 1;

% 遍历每个股票符号并获取历史数据
for i = 1:length(tickers)
    % 更新等待条显示当前股票
    waitbar((i-1)/length(tickers), h, ...
        sprintf('正在获取 %s 的历史股票数据 (%0.2f%%)', ...
        tickers{i}, (i-1)*100/length(tickers)))

    % 创建用于检索数据的 URL 字符串
    url = sprintf(['https://query1.finance.yahoo.com/v7/finance/download/', ...
        '%s?period1=%d&period2=%d&interval=1%s&events=%s'], ...
        tickers{i}, startDate, endDate, freq, event);

    % 从 Yahoo Finance 获取数据
    [temp, status] = urlread(url,'post',{'matlabstockdata@yahoo.com', 'historical stocks'});
        
    % 如果成功下载数据，则处理它。否则，忽略此股票符号
    if status
        % 解析历史数据
        [date, op, high, low, cl, adj_close, volume] = ...
            strread(temp(43:end),'%s%s%s%s%s%s%s','delimiter',',');
            
        % 存储数据
        stocks(idx).Ticker = tickers{i};
        stocks(idx).Date = date;
        stocks(idx).Open = str2double(op);
        stocks(idx).High = str2double(high);
        stocks(idx).Low = str2double(low);
        stocks(idx).Close = str2double(cl);
        stocks(idx).AdjClose = str2double(adj_close);
        stocks(idx).Volume = str2double(volume);
        
        idx = idx + 1;  % 增加股票索引
    end
    
    % 更新等待条
    waitbar(i/length(tickers), h)
end

close(h)  % 关闭等待条
save('stocks_restt.mat', 'stocks')
