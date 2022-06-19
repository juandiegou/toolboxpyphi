function mv_qqplot(data)

% -----------------------------------------------------------------------
% This function produces a QQ-plot to assess multivariate normality. If
% your data are multivariate normal, then there should be a linear
% relationship between the ordered Mahalanobis distances of the data from
% the mean vector and their corresponding chi-square quantiles.
% Graphically, this function should produce a series of blue points roughly
% along the red line if your data are normal.
% -----------------------------------------------------------------------

% find where this code was modified from 

% just take the amount of data that won't make Matlab crap out
if size(data,2)<2000 && size(data,1)<50
    X=data';
elseif size(data,2)<2000 && size(data,1)>=50
    X=data(1:50,:)';
elseif size(data,2)>=2000 && size(data,1)<50
    X=data(:,1:2000)';
else
    X=data(1:50,1:2000)';
end



[n,p] = size(X);

difT = [];
for	j = 1:p
   eval(['difT=[difT,(X(:,j)-mean(X(:,j)))];']);
end
S=cov(X);
D = difT*inv(S)*difT'; 
[d,t] = sort(diag(D));   %squared Mahalanobis distances
r = tiedrank(d); 
chi2q=chi2inv((r-0.5)./n,p);

plot(d,chi2q,'*');
hold on
e=round(max([max(d) max(chi2q)]));
plot(1:e,1:e,'r')
xlabel('Mahalanobis distance quantile');
ylabel('Chi-square quantile')

xlim([min(d) round(e*1.1)])
