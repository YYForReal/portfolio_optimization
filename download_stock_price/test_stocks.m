% 检查一下数据
% 假设 stocks 是结构体数组

stocks = load("clean_FTSE100_stocks.mat").new_stocks;
fprintf("length(stocks) %d\n",length(stocks))
for i = 1:min(10, length(stocks))
    ticker = stocks(i).Ticker;
    dates = stocks(i).Date;
    fprintf('Ticker: %s, Date Length: %d\n', ticker, length(dates));
end
