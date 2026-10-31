function test_qsvds_backends
%TEST_QSVDS_BACKENDS Cross-backend accuracy, structure, status and options.
rng(321,'twister'); count=0;
for backend={'svds','krylov','subspace'}
    for dims={[9,7],[7,9],[12,12]}
        d=dims{1}; m=d(1); n=d(2); p=min(m,n);
        matrices={quaternion_to_complex(randn(m,n)+1i*randn(m,n),randn(m,n)+1i*randn(m,n))};
        deficient=zeros(m,n); deficient(1:2,1:2)=diag([8,3]);
        matrices{end+1}=quaternion_to_complex(deficient,zeros(m,n));
        if m==n
            matrices{end+1}=quaternion_to_complex(diag([8,8,8,8,3,2,1,zeros(1,5)]),zeros(12));
            matrices{end+1}=quaternion_to_complex(eye(12),zeros(12));
        end
        for c=1:numel(matrices)
            A=matrices{c}; s=svd(A); s=s(1:2:2*p);
            for k=unique([1,3,p-1,p])
                state=rng;
                [U,E,V,info]=QSVDS(A,m,n,k,struct('backend',backend{1}));
                assert(isequal(state,rng),'RNG was not restored.');
                assert(info.flag==0 && info.converged,'Failure: %s, %dx%d, case %d, k=%d.',backend{1},m,n,c,k);
                assert(isequal(size(U),[2*m,2*k]) && isequal(size(V),[2*n,2*k]));
                assert(norm(diag(E(1:k,1:k))-s(1:k))<1e-8*max(1,s(1)));
                assert(norm(A*V-U*E,'fro')<1e-8*max(1,norm(A,'fro')));
                assert(norm(A'*U-V*E','fro')<1e-8*max(1,norm(A,'fro')));
                assert(max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro'))<1e-8);
                assert(norm(U-quaternion_to_complex(U(1:m,1:k),U(1:m,k+1:2*k)),'fro')<1e-11);
                assert(norm(V-quaternion_to_complex(V(1:n,1:k),V(1:n,k+1:2*k)),'fro')<1e-11);
                if strcmp(backend{1},'krylov'), assert(~info.used_builtin_svds); end
                if k==p, assert(info.full_decomposition && ~info.used_builtin_svds); end
                count=count+1;
            end
        end
    end
end
A=quaternion_to_complex(eye(4),zeros(4));
for opts={struct('backend','bad'),struct('maxspace',0),struct('tol',NaN),struct('maxit',0),struct('seed',-1),struct('extra',1)}
    failed=false; try, QSVDS(A,4,4,1,opts{1}); catch, failed=true; end
    assert(failed,'Invalid option accepted.');
end
% A flat spectrum with distinct close singular values must now converge.
A=quaternion_to_complex(diag(linspace(1,.8,120)),zeros(120));
[U,E,V,info]=QSVDS(A,120,120,3,struct('backend','krylov'));
assert(info.flag==0 && ~info.used_builtin_svds && info.iterations<1000);
assert(norm(A*V-U*E,'fro')<1e-8 && norm(A'*U-V*E','fro')<1e-8);
fprintf('All %d backend accuracy checks, invalid-option checks and pure-Krylov flat-spectrum check passed.\n',count);
end
