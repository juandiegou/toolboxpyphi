function Cov_X_Y = cov_cond(Cov_X,Cov_XY,Cov_Y)

%-----------------------------------------------------------------------
% This function calculates the partial covariance of the present of your
% system given its past. Its inputs are the covariance of the present
% (Cov_X), the cross-covariance between the present and the past (Cov_XY),
% and the covariance of the past (Cov_Y).
%-----------------------------------------------------------------------

Cov_X_Y = Cov_X - Cov_XY/Cov_Y*Cov_XY';
