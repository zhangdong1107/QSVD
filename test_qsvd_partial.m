function test_qsvd_partial
%TEST_QSVD_PARTIAL Compare top-K values/factors with a full complex SVD.
rng(123,'twister');
checks = 0;
shapes = {[12,9],[9,12],[8,8]};
for t = 1:numel(shapes)
    dims = shapes{t}; m = dims(1); n = dims(2); p = min(m,n);
    AC = quaternion_to_complex(randn(m,n)+1i*randn(m,n),randn(m,n)+1i*randn(m,n));
    reference = svd(AC); reference = reference(1:2:2*p);
    for k = unique([1,3,p-1,p])
        [s,info] = qsvd_values(AC,m,n,k);
        assert(isequal(size(s),[k,1]) && isreal(s) && all(diff(s)<=0));
        assert(norm(s-reference(1:k)) < 1e-9*norm(reference));
        if k < p, assert(info.flag == 0 && ~info.full_decomposition && info.largest_projected_svd < p); end
        s2 = qsvd_complex(AC,m,n,[],k);
        assert(norm(s2-s) < 1e-9*norm(reference));
        [U,E,V] = qsvd_complex(AC,m,n,[],k);
        assert(isequal(size(U),[2*m,2*k]) && isequal(size(E),[2*k,2*k]) && isequal(size(V),[2*n,2*k]));
        assert(norm(U'*U-eye(2*k),'fro') < 1e-9);
        assert(norm(V'*V-eye(2*k),'fro') < 1e-9);
        assert(norm(AC*V-U*E,'fro') < 1e-9*norm(AC,'fro'));
        assert(norm(AC'*U-V*E','fro') < 1e-9*norm(AC,'fro'));
        % Optimal rank-K quaternion error = doubled complex tail energy.
        assert(abs(norm(AC-U*E*V','fro')-sqrt(2*sum(reference(k+1:end).^2))) < 1e-9*norm(AC,'fro'));
        assert(norm(diag(E(1:k,1:k))-s) < 1e-9*norm(reference));
        assert(norm(U-quaternion_to_complex(U(1:m,1:k),U(1:m,k+1:2*k)),'fro') < 1e-12);
        assert(norm(V-quaternion_to_complex(V(1:n,1:k),V(1:n,k+1:2*k)),'fro') < 1e-12);
        checks = checks+1;
    end
    assert(abs(qsvd_values(AC,m,n)-reference(1)) < 1e-9*norm(reference));
    % Four-argument, single-output behavior is still the left factor.
    Uold = qsvd_complex(AC,m,n,1e-12);
    assert(isequal(size(Uold),[2*m,2*p]));
end
% Repetition across the requested cutoff, zeros, scalar, and scale extremes.
cases = {diag([7,7,7,3,1,0,0,0]), diag([8,2,0,0,0,0,0,0]), zeros(5,3), zeros(3,5), 2+3i};
for c = 1:numel(cases)
    B = cases{c}; [m,n] = size(B); p = min(m,n);
    AC = quaternion_to_complex(B,zeros(m,n));
    reference = svd(AC); reference = reference(1:2:2*p);
    for k = unique([1,min(2,p),p])
        s = qsvd_values(AC,m,n,k);
        assert(norm(s-reference(1:k)) < 1e-9*max(1,norm(reference)));
        [U,E,V] = qsvd_complex(AC,m,n,[],k);
        assert(isequal(size(U),[2*m,2*k]) && isequal(size(V),[2*n,2*k]));
        assert(max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro')) < 1e-9);
        assert(norm(AC*V-U*E,'fro') < 1e-9*max(1,norm(AC,'fro')));
        checks = checks+1;
    end
end
for scale = [1e-150,1e150]
    AC = quaternion_to_complex(scale*diag(9:-1:1),zeros(9));
    assert(abs(qsvd_values(AC,9,9,1)/scale-9) < 1e-9);
    checks = checks+1;
end
AC = quaternion_to_complex(eye(3),zeros(3));
for k = {0,-1,1.5,4,NaN,Inf,[1,2],1+1i}
    assert_fails(@() qsvd_values(AC,3,3,k{1}));
    assert_fails(@() qsvd_complex(AC,3,3,[],k{1}));
    checks = checks+1;
end
assert_fails(@() qsvd_values(eye(5),3,3,1));
bad = AC; bad(1,4) = 2;
assert_fails(@() qsvd_values(bad,3,3,1));
bad = AC; bad(1) = Inf;
assert_fails(@() qsvd_values(bad,3,3,1));
assert_fails(@() qsvd_complex(AC,3,3,-1,1));
fprintf('All %d top-K QSVD checks passed, plus input validation.\n',checks);
end

function assert_fails(f)
failed = false;
try, f(); catch, failed = true; end
assert(failed,'Invalid input was not rejected.');
end
