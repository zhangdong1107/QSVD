function [UC, EC, VC] = QSVD(AC,m,n,tol)
%QSVD Full quaternion SVD using the author's original complex algorithm.
%   [UC,EC,VC] = QSVD(AC,M,N,TOL) returns complex-representation factors.
%   S = QSVD(AC,M,N,TOL) returns all min(M,N) quaternion singular values.
%   TOL defaults to 1e-12. See QSVDS for an iterative partial decomposition.
if nargin < 4, tol = []; end
[UC,EC,VC] = qsvd_complex(AC,m,n,tol);
if nargout <= 1
    p = min(m,n);
    UC = diag(EC(1:p,1:p));
end
end
