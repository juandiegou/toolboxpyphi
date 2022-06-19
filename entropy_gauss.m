function H = entropy_gauss(Cov_X)

%-----------------------------------------------------------------------
% This function calculates the differential entropy of a multivariate
% normal signal from its covariance matrix.
%-----------------------------------------------------------------------


n = size(Cov_X,1);

if isinf(det(Cov_X))
    % this is a little trick for getting the log of the determinant of
    % really large matrices, which can sometimes incorrectly give you
    % infinity
    L = chol(Cov_X);
    logdet_Cov_X = 2*sum(log(diag(L)));   
    H = 1/2*logdet_Cov_X + 1/2*n*log(2*pi*exp(1));
else
    H = 1/2*log(det(Cov_X)) + 1/2*n*log(2*pi*exp(1));
end