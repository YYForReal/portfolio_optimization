function get_weekly_price

% load('stocks_weekly_FTSE.mat');
load("stoskc.mat")

for i=1:size(new_stocks,2)
    wk_price(:,i)=new_stocks(i).AdjClose;
    % wk_return(:,i) = price2ret(new_stocks(i).AdjClose,[],'periodic');
    % wk_return(:,i) = price2ret(new_stocks(i).AdjClose);
    wk_return(:, i) = diff(new_stocks(i).AdjClose) ./ new_stocks(i).AdjClose(1:end-1);

    mean_return(i,:)=mean(wk_return(:,i));
    stdDev_return(i,:)=std(wk_return(:,i));
end
Correlation=corrcoef(wk_return);

save('correlation_matrix', 'Correlation','-ASCII')
save('mean_return', 'mean_return','-ASCII')
save('stdDev_return', 'stdDev_return','-ASCII')
save('wk_price', 'wk_price','-ASCII')
save('wk_return', 'wk_return','-ASCII')

date=new_stocks(1).Date;
dn=datenum(date,'yyyy-mm-dd');
save('date', 'dn','-ASCII');

end