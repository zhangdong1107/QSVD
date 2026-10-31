function test_qsvds_real_reference(backend)
%TEST_QSVDS_REAL_REFERENCE Independent 4-by-4 real representation checks.
% This reference uses neither the full QSVD nor the doubled complex SVD.
if nargin<1, backend='svds'; end
rng(2026,'twister');
checks = 0; worst_value = 0; worst_residual = 0; worst_orth = 0;
for dims = {[7,4],[4,7],[9,9],[15,11]}
    d = dims{1}; m = d(1); n = d(2); p = min(m,n);
    A0 = randn(m,n); A1 = randn(m,n); A2 = randn(m,n); A3 = randn(m,n);
    AC = quaternion_to_complex(A0+1i*A1,A2+1i*A3);
    AR = real_representation(A0,A1,A2,A3);
    reference = svd(AR); reference = reference(1:4:4*p);
    for seed = [0,7,31]
        for k = unique([1,min(3,p-1),p-1])
            [U,E,V,info] = QSVDS(AC,m,n,k,struct('seed',seed,'backend',backend));
            assert(info.converged && ~info.full_decomposition);
            values = diag(E(1:k,1:k));
            UR = from_complex(U,m,k); VR = from_complex(V,n,k);
            ER = real_representation(diag(values),zeros(k),zeros(k),zeros(k));
            value_error = norm(values-reference(1:k))/norm(reference(1:k));
            residual = max(norm(AR*VR-UR*ER,'fro'),norm(AR'*UR-VR*ER,'fro'))/norm(AR,'fro');
            orth = max(norm(UR'*UR-eye(4*k),'fro'),norm(VR'*VR-eye(4*k),'fro'));
            assert(value_error < 1e-9 && residual < 1e-9 && orth < 1e-9);
            worst_value = max(worst_value,value_error);
            worst_residual = max(worst_residual,residual);
            worst_orth = max(worst_orth,orth);
            checks = checks+1;
        end
    end
end
% Analytic rank-one singular value = ||x||_H ||y||_H.
x = randn(7,4); y = randn(5,4);
XC = quaternion_to_complex(x(:,1)+1i*x(:,2),x(:,3)+1i*x(:,4));
YC = quaternion_to_complex(y(:,1)+1i*y(:,2),y(:,3)+1i*y(:,4));
[~,E,~,info] = QSVDS(XC*YC',7,5,1,struct('backend',backend));
expected = norm(x,'fro')*norm(y,'fro');
assert(info.converged && abs(E(1,1)-expected)/expected < 1e-12);
fprintf('All %d independent real-reference QSVDS checks and analytic rank-one check passed.\n',checks);
fprintf('Worst relative value error %.3e; real singular-equation residual %.3e; real orthogonality %.3e.\n',worst_value,worst_residual,worst_orth);
end

function R = real_representation(A0,A1,A2,A3)
R = [A0,-A1,-A2,-A3; A1,A0,-A3,A2; A2,A3,A0,-A1; A3,-A2,A1,A0];
end

function R = from_complex(C,m,n)
B1 = C(1:m,1:n); B2 = C(1:m,n+1:2*n);
R = real_representation(real(B1),imag(B1),real(B2),imag(B2));
end
