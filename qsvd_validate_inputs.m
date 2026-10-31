function p = qsvd_validate_inputs(AC, m, n, k)
%QSVD_VALIDATE_INPUTS Shared validation for the QSVD entry points.
validateattributes(m, {'double'}, {'real','scalar','finite','integer','positive'}, mfilename, 'm');
validateattributes(n, {'double'}, {'real','scalar','finite','integer','positive'}, mfilename, 'n');
validateattributes(AC, {'double'}, {'2d','finite','nonsparse','size',[2*m,2*n]}, mfilename, 'AC');
p = min(m,n);
if nargin == 4
    validateattributes(k, {'numeric'}, {'real','scalar','finite','integer','>=',1,'<=',p}, mfilename, 'k');
end
B1 = AC(1:m,1:n); B2 = AC(1:m,n+1:2*n);
if norm(AC-quaternion_to_complex(B1,B2),'fro') > 1e-12*max(1,norm(AC,'fro'))
    error('qsvd:InvalidRepresentation', 'AC must have blocks [B1,B2;-conj(B2),conj(B1)].');
end
end
