function stocks = hist_stock_data(start_date, end_date, varargin)
% �������ڻ�ȡ��ʷ��Ʊ����

stocks = struct([]);        % ��ʼ�����ݽṹ

%% ��������

% ��ʽ����ʼ�ͽ�������Ϊ Posix ʱ��
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

% ȷ���û��Ƿ�ָ����Ƶ��
temp = find(strcmp(varargin,'frequency') == 1);
 if isempty(temp)
    freq = 'd';  % ���δ������Ĭ��Ϊÿ������
else
    freq = varargin{temp+1};
    varargin(temp:temp+1) = [];
end
clear temp

% ȷ���û��Ƿ�ָ�����¼�����
temp = find(strcmp(varargin,'type') == 1);
if isempty(temp)
    event = 'history';  % Ĭ��Ϊ��ʷ�۸�
else
    event = varargin{temp+1};
    varargin(temp:temp+1) = [];
end
clear temp

% ��� varargin �ĵ�һ��Ԫ�ر�����һ����Ԫ���飬��ٶ�����һ����Ʊ���ŵĵ�Ԫ����
if iscell(varargin{1})
    tickers = varargin{1};
% ���򣬼�����Ƿ��� .txt �ļ�
elseif ~isempty(strfind(varargin{1},'.txt'))
    fid = fopen(varargin{1}, 'r');
    tickers = textscan(fid, '%s'); tickers = tickers{:};
    fclose(fid);
% ���򣬼ٶ���Ҫô�ǵ�����Ʊ���ţ�Ҫô�ǹ�Ʊ�����б�
else
    tickers = varargin;
end

%% ��ȡ��ʷ����

h = waitbar(0, '��ȴ�...');  % �����ȴ���
idx = 1;

% ����ÿ����Ʊ���Ų���ȡ��ʷ����
for i = 1:length(tickers)
    % ���µȴ�����ʾ��ǰ��Ʊ
    waitbar((i-1)/length(tickers), h, ...
        sprintf('���ڻ�ȡ %s ����ʷ��Ʊ���� (%0.2f%%)', ...
        tickers{i}, (i-1)*100/length(tickers)))

    % �������ڼ������ݵ� URL �ַ���
    url = sprintf(['https://query1.finance.yahoo.com/v7/finance/download/', ...
        '%s?period1=%d&period2=%d&interval=1%s&events=%s'], ...
        tickers{i}, startDate, endDate, freq, event);

    % �� Yahoo Finance ��ȡ����
    [temp, status] = urlread(url,'post',{'matlabstockdata@yahoo.com', 'historical stocks'});
        
    % ����ɹ��������ݣ������������򣬺��Դ˹�Ʊ����
    if status
        % ������ʷ����
        [date, op, high, low, cl, adj_close, volume] = ...
            strread(temp(43:end),'%s%s%s%s%s%s%s','delimiter',',');
            
        % �洢����
        stocks(idx).Ticker = tickers{i};
        stocks(idx).Date = date;
        stocks(idx).Open = str2double(op);
        stocks(idx).High = str2double(high);
        stocks(idx).Low = str2double(low);
        stocks(idx).Close = str2double(cl);
        stocks(idx).AdjClose = str2double(adj_close);
        stocks(idx).Volume = str2double(volume);
        
        idx = idx + 1;  % ���ӹ�Ʊ����
    end
    
    % ���µȴ���
    waitbar(i/length(tickers), h)
end

close(h)  % �رյȴ���
save('stocks_restt.mat', 'stocks')
