function [s, info] = qsvd_values(AC,m,n,k)
%QSVD_VALUES Compatibility helper for the iterative QSVDS values.
%   Delegates to QSVDS with the default MATLAB SVDS backend.
%   K defaults to 1. Inspect INFO.flag and INFO.relative_residuals.
if nargin < 4 || isempty(k), k = 1; end
[~,E,~,info] = QSVDS(AC,m,n,k);
s = diag(E(1:k,1:k));
end
