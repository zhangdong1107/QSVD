function test_qsvds
%TEST_QSVDS Numerical and execution-path checks for the iterative solver.
rng(87,'twister');
checks = 0;
spectra = {linspace(15,1,20), [9,9,9,9,5,4,3,2,ones(1,12)], ...
    [10,9.9999,9.9998,9.9997,8:-1:1,zeros(1,8)], ...
    [8,4,2,zeros(1,17)], ones(1,20)};
for c = 1:numel(spectra)
    s = spectra{c}; n = numel(s); m = n+4;
    [L,~,~] = QSVD(quaternion_to_complex(randn(m)+1i*randn(m),randn(m)+1i*randn(m)),m,m);
    [R,~,~] = QSVD(quaternion_to_complex(randn(n)+1i*randn(n),randn(n)+1i*randn(n)),n,n);
    D = zeros(m,n); D(1:n,1:n) = diag(s);
    AC = L*quaternion_to_complex(D,zeros(m,n))*R';
    for k = [1,3,5]
        state = rng;
        [U,E,V,info] = QSVDS(AC,m,n,k);
        assert(isequal(state,rng),'QSVDS changed the caller random state.');
        assert(info.flag == 0 && info.converged);
        assert(~info.full_decomposition && info.largest_projected_svd < min(m,n));
        assert(norm(diag(E(1:k,1:k))-s(1:k).') < 2e-8*max(s));
        assert(max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro')) < 1e-9);
        residual = max(norm(AC*V-U*E,'fro'),norm(AC'*U-V*E','fro'))/norm(AC,'fro');
        assert(residual < 1e-9,'Singular-vector residual too large.');
        assert(max(info.relative_residuals) <= 1e-10);
        checks = checks+1;
    end
end
% Explicit nonconvergence must remain iterative and visible to the caller.
AC = quaternion_to_complex(randn(30,25)+1i*randn(30,25),randn(30,25)+1i*randn(30,25));
saved_warning = warning('off','qsvds:NoConvergence');
restore_warning = onCleanup(@() warning(saved_warning)); %#ok<NASGU>
lastwarn('');
[~,~,~,info] = QSVDS(AC,30,25,2,struct('maxit',1,'tol',1e-14,'backend','krylov'));
assert(info.flag == 1 && ~info.converged && ~info.full_decomposition && info.iterations == 1);
assert(info.largest_projected_svd < 25);
% Also check a genuinely larger low-rank matrix and a wide zero matrix.
X = quaternion_to_complex(randn(160,8)+1i*randn(160,8),randn(160,8)+1i*randn(160,8));
Y = quaternion_to_complex(randn(120,8)+1i*randn(120,8),randn(120,8)+1i*randn(120,8));
ACbig = X*Y';
[U,E,V,info] = QSVDS(ACbig,160,120,3);
reference = svd(ACbig); reference = reference(1:2:6);
assert(info.converged && info.largest_projected_svd < 120);
assert(norm(diag(E(1:3,1:3))-reference) < 1e-8*norm(reference));
assert(norm(ACbig*V-U*E,'fro')/norm(ACbig,'fro') < 1e-9);
[U,E,V,info] = QSVDS(zeros(6,10),3,5,2);
assert(isequal(size(U),[6,4]) && isequal(size(V),[10,4]) && ~any(E(:)) && info.flag == 0);
% Runtime guard: forbidden iterative built-ins must never be called.
guard_dir = tempname; mkdir(guard_dir);
for name = {'svds','eigs'}
    file = fopen(fullfile(guard_dir,[name{1},'.m']),'w');
    fprintf(file,'function varargout = %s(varargin)\nerror(''qsvds:ForbiddenBuiltin'',''Forbidden built-in called.'');\nend\n',name{1});
    fclose(file);
end
% MATLAB can bind built-ins before path shadowing. Instrument copies of the
% two production files at the actual SVD call sites instead. Only the names
% and the SVD call target change; compare their outputs with the originals.
core = fileread(which('qsvd_complex'));
core = strrep(core,'qsvd_complex','qsvds_audited_core');
core = regexprep(core,'\<svd\s*\(','qsvds_checked_svd(');
partial = fileread(which('qsvds_krylov'));
partial = strrep(partial,'qsvds_krylov','qsvds_audited_partial');
partial = strrep(partial,'qsvd_complex','qsvds_audited_core');
file = fopen(fullfile(guard_dir,'qsvds_audited_core.m'),'w'); fprintf(file,'%s',core); fclose(file);
file = fopen(fullfile(guard_dir,'qsvds_audited_partial.m'),'w'); fprintf(file,'%s',partial); fclose(file);
file = fopen(fullfile(guard_dir,'qsvds_checked_svd.m'),'w');
fprintf(file,['function varargout = qsvds_checked_svd(varargin)\n' ...
    'global QSVDS_TEST_FULL_DIM QSVDS_TEST_MAX_DIM\n' ...
    'd=min(size(varargin{1})); QSVDS_TEST_MAX_DIM=max(QSVDS_TEST_MAX_DIM,d);\n' ...
    'if d>=QSVDS_TEST_FULL_DIM, error(''qsvds:ForbiddenFullSVD'',''Full-input SVD called.''); end\n' ...
    '[varargout{1:nargout}]=builtin(''svd'',varargin{:});\nend\n']);
fclose(file);
global QSVDS_TEST_FULL_DIM QSVDS_TEST_MAX_DIM
QSVDS_TEST_FULL_DIM = 2*25; QSVDS_TEST_MAX_DIM = 0;
addpath(guard_dir,'-begin');
guard_cleanup = onCleanup(@() cleanup_guard(guard_dir)); %#ok<NASGU>
rehash;
[U0,E0,V0,info0] = QSVDS(AC,30,25,2,struct('backend','krylov'));
[U1,E1,V1,info] = qsvds_audited_partial(AC,30,25,2);
assert(info.converged && info0.converged);
assert(norm(U0-U1,'fro')+norm(E0-E1,'fro')+norm(V0-V1,'fro') < 1e-10);
assert(QSVDS_TEST_MAX_DIM > 0 && QSVDS_TEST_MAX_DIM <= 2*info.largest_projected_svd, ...
    'Observed SVD dimension %g; expected <= %g.',QSVDS_TEST_MAX_DIM,2*info.largest_projected_svd);
fprintf('Largest intercepted complex SVD dimension %d (full input minimum dimension 50).\n',QSVDS_TEST_MAX_DIM);
fprintf('All %d structured QSVDS spectrum checks passed; large/zero/nonconvergence/RNG/forbidden-builtin checks passed.\n',checks);
end

function cleanup_guard(directory)
rmpath(directory);
clear qsvds_audited_partial qsvds_audited_core qsvds_checked_svd
for name = {'svds.m','eigs.m','qsvds_audited_core.m','qsvds_audited_partial.m','qsvds_checked_svd.m'}
    delete(fullfile(directory,name{1}));
end
rmdir(directory);
clear global QSVDS_TEST_FULL_DIM QSVDS_TEST_MAX_DIM
end
