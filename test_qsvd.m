function test_qsvd
%TEST_QSVD Check reconstruction, orthogonality, structure and singular values.
%   Run after extraction: test_qsvd

rng(42, 'twister');
cases = cell(0,3);
for shape = {[4,4], [5,3], [3,5]}
    dims = shape{1}; m = dims(1); n = dims(2);
    B1 = randn(m,n) + 1i*randn(m,n);
    B2 = randn(m,n) + 1i*randn(m,n);
    cases(end+1,:) = {sprintf('random_%dx%d',m,n), B1, B2};
end
cases(end+1,:) = {'repeated_values', diag([4,4,1,1]), zeros(4)};
cases(end+1,:) = {'rank_deficient', diag([4,1,0,0]), zeros(4)};
cases(end+1,:) = {'zero_tall', zeros(4,2), zeros(4,2)};
cases(end+1,:) = {'zero_wide', zeros(2,4), zeros(2,4)};
cases(end+1,:) = {'scalar', 2+3i, 4-1i};

for c = 1:size(cases,1)
    name = cases{c,1}; B1 = cases{c,2}; B2 = cases{c,3};
    [m,n] = size(B1); p = min(m,n);
    AC = quaternion_to_complex(B1,B2);
    [UC,EC,VC] = qsvd_complex(AC,m,n);
    reconstruction = norm(AC - UC*EC*VC','fro') / max(1,norm(AC,'fro'));
    orthogonality = max(norm(UC'*UC-eye(size(UC,2)),'fro'), ...
                        norm(VC'*VC-eye(size(VC,2)),'fro'));
    values = svd(AC);
    actual = diag(EC(1:p,1:p));
    value_error = norm(actual-values(1:2:2*p)) / max(1,norm(values));
    pu = size(UC,2)/2; pv = size(VC,2)/2;
    u_structure = norm(UC-quaternion_to_complex(UC(1:m,1:pu),UC(1:m,pu+1:2*pu)),'fro');
    v_structure = norm(VC-quaternion_to_complex(VC(1:n,1:pv),VC(1:n,pv+1:2*pv)),'fro');
    if norm(AC,'fro') == 0
        assert(isequal(size(UC),[2*m,2*m]) && isequal(size(EC),[2*m,2*n]) && isequal(size(VC),[2*n,2*n]));
    else
        assert(isequal(size(UC),[2*m,2*p]) && isequal(size(EC),[2*p,2*p]) && isequal(size(VC),[2*n,2*p]));
    end
    assert(max([reconstruction,orthogonality,value_error,u_structure,v_structure]) < 1e-10, ...
        'QSVD check failed for %s.', name);
    fprintf('PASS %-20s residual %.3e  orthogonality %.3e\n',name,reconstruction,orthogonality);
end
fprintf('All %d QSVD checks passed.\n',size(cases,1));
end
