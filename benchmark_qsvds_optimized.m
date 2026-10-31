function benchmark_qsvds_optimized(mode,output_dir)
%BENCHMARK_QSVDS Wall time and accuracy, known-spectrum quaternion matrices.
%   Run from the homepage checkout; the installed solver is not modified.
if nargin < 1, mode = 'full'; end
script_dir = fileparts(mfilename('fullpath'));
if exist(fullfile(script_dir,'solver','QSVDS.m'),'file')
    solver_dir = fullfile(script_dir,'solver');
    if nargin < 2, output_dir = script_dir; end
elseif exist(fullfile(script_dir,'QSVD.m'),'file')
    solver_dir = script_dir;
    if nargin < 2, output_dir = fullfile(script_dir,'benchmark-results'); end
else
    solver_dir = fullfile(fileparts(script_dir),'matlab','QSVD-AML-MATLAB');
    if nargin < 2, output_dir = fullfile(fileparts(fileparts(script_dir)),'qsvds-optimization-20261031'); end
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
addpath(solver_dir);
tol = 1e-10; maxit = 1000;
if strcmp(mode,'pilot')
    dims = [256,192]; repeats = 1;
else
    dims = [512,384;1024,768;1536,1024]; repeats = 3;
end
ks = [1,5]; spectra = {'low_rank_dominant','power_decay','flat'};
metadata = struct('matlab',version,'release',version('-release'), ...
    'computer',computer,'threads',maxNumCompThreads,'mode',mode, ...
    'tolerance',tol,'maxit',maxit,'repeats',repeats,'dimensions',dims, ...
    'k',ks,'time_unit','seconds','reference','Known analytic spectrum after quaternion Householder rotations', ...
    'timing','Matrix generation, reference construction and accuracy checks excluded; factors included', ...
    'methods',{{'QSVD','QSVDS_krylov','QSVDS','MATLAB_svds'}},'default_QSVDS','SVDS + paired recovery and thin refinement; internal tolerance=1e-11', ...
    'native_svds','Requests 2*k complex triplets; no quaternion factor conversion included');
write_json(fullfile(output_dir,[mode,'-metadata.json']),metadata);
fprintf('MATLAB %s; threads %d; %s, repeats %d\n',version,metadata.threads,mode,repeats);
% Warm all three paths before measuring; this is not a timed result.
[warm,~] = make_matrix(64,48,'power_decay',30);
[~,~,~] = QSVD(warm,64,48);
[~,~,~,~] = QSVDS(warm,64,48,2);
[~,~,~,~] = svds(warm,4,'largest',struct('tol',tol,'maxit',maxit));
prototype = make_row('',0,0,0,'',0,0,0,0,0, ...
    struct('value_error',0,'residual',0,'orthogonality',0,'pair_error',0));
rows = repmat(prototype,0,1); full_times = zeros(repeats,1);
saved_warning = warning('off','qsvds:NoConvergence');
warning_cleanup = onCleanup(@() warning(saved_warning)); %#ok<NASGU>
for shape = 1:size(dims,1)
    m = dims(shape,1); n = dims(shape,2); p = min(m,n);
    for c = 1:numel(spectra)
        family = spectra{c};
        [AC,truth] = make_matrix(m,n,family,100*shape+c);
        fprintf('\nCASE %dx%d quaternion (%dx%d complex), %s\n',m,n,2*m,2*n,family);
        full_metrics = cell(numel(ks),repeats);
        % Full QSVD is independent of k; measure it once per repeat and reuse
        % that cost for both k comparisons. Partial method order alternates.
        for rep = 1:repeats
            rng(900+rep,'twister');
            timer = tic; [U,E,V] = QSVD(AC,m,n); full_times(rep) = toc(timer);
            for a = 1:numel(ks)
                k = ks(a); columns = [1:k,p+(1:k)];
                full_metrics{a,rep} = measure_accuracy(AC,U(:,columns),E(columns,columns),V(:,columns),truth(1:k));
            end
            fprintf('  full QSVD rep %d: %.4f s\n',rep,full_times(rep)); drawnow;
        end
        clear U E V
        for a = 1:numel(ks)
            k = ks(a);
            for rep = 1:repeats
                rows(end+1) = make_row(family,m,n,k,'QSVD',rep,full_times(rep),0,0,p,full_metrics{a,rep}); %#ok<AGROW>
                order = {'QSVDS_krylov','QSVDS','MATLAB_svds'};
                if mod(rep,2)==0, order = fliplr(order); end
                for b = 1:numel(order)
                    method = order{b}; rng(900+rep,'twister');
                    if ~strcmp(method,'MATLAB_svds')
                        timer = tic;
                        callopts=struct('tol',tol,'maxit',maxit);
                        if strcmp(method,'QSVDS_krylov'), callopts.backend='krylov'; end
                        [U,E,V,info] = QSVDS(AC,m,n,k,callopts);
                        elapsed = toc(timer); flag = info.flag; iterations = info.iterations; ell = info.subspace_dimension;
                    else
                        timer = tic;
                        [U,E,V,flag] = svds(AC,2*k,'largest',struct('tol',tol,'maxit',maxit));
                        elapsed = toc(timer); iterations = NaN; ell = NaN;
                    end
                    metrics = measure_accuracy(AC,U,E,V,truth(1:k));
                    metrics.actual_backend=method; metrics.used_builtin_svds=strcmp(method,'MATLAB_svds'); metrics.fallback_reason='';
                    if ~strcmp(method,'MATLAB_svds'), metrics.actual_backend=info.method; metrics.used_builtin_svds=info.used_builtin_svds; metrics.fallback_reason=info.fallback_reason; end
                    rows(end+1) = make_row(family,m,n,k,method,rep,elapsed,flag,iterations,ell,metrics); %#ok<AGROW>
                    fprintf('  k=%d %-12s rep %d: %.4f s flag %d iter %g value %.2e residual %.2e orth %.2e\n', ...
                        k,method,rep,elapsed,flag,iterations,metrics.value_error,metrics.residual,metrics.orthogonality);
                    write_json(fullfile(output_dir,[mode,'-raw.json']),rows); drawnow;
                    clear U E V
                end
            end
        end
        clear AC
    end
end
data = struct2table(rows);
writetable(data,fullfile(output_dir,[mode,'-raw.csv']));
save(fullfile(output_dir,[mode,'-results.mat']),'metadata','rows');
fprintf('\nBENCHMARK COMPLETE: %s\n',output_dir);
end

function [A,s] = make_matrix(m,n,family,seed)
rng(seed,'twister'); p = min(m,n);
switch family
    case 'low_rank_dominant'
        s = [10,8,6,4,3,2,1.5,1,1e-3*linspace(.99,.01,p-8)].';
    case 'power_decay'
        s = (1./(1:p).^2).';
    case 'flat'
        s = linspace(1,.8,p).';
end
D = zeros(m,n); D(1:p,1:p) = diag(s);
A = quaternion_to_complex(D,zeros(m,n));
% Independent quaternion unitary Householder rotations preserve the known
% spectrum. They make the dense matrix have nonzero quaternion components.
% Each H=I-2*W*W', W=Chi(w), W'*W=I_2; no reference QSVD is needed.
for j = 1:16
    x = randn(m,4); x = x/norm(x,'fro');
    W = quaternion_to_complex(x(:,1)+1i*x(:,2),x(:,3)+1i*x(:,4));
    A = A-2*W*(W'*A);
    y = randn(n,4); y = y/norm(y,'fro');
    Z = quaternion_to_complex(y(:,1)+1i*y(:,2),y(:,3)+1i*y(:,4));
    A = A-2*(A*Z)*Z';
end
end

function metrics = measure_accuracy(A,U,E,V,truth)
k = numel(truth);
actual = sort(real(diag(E)),'descend');
expected = repelem(truth,2);
metrics.value_error = norm(actual-expected)/norm(expected);
metrics.pair_error = max(abs(actual(1:2:end)-actual(2:2:end)))/truth(1);
R1 = A*V-U*E; R2 = A'*U-V*E';
metrics.residual = max(sqrt(sum(abs(R1).^2,1)+sum(abs(R2).^2,1)))/truth(1);
metrics.orthogonality = max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro'));
end

function row = make_row(family,m,n,k,method,rep,seconds,flag,iterations,ell,metrics)
row = struct('family',family,'m',m,'n',n,'k',k,'method',method,'repeat',rep, ...
    'seconds',seconds,'flag',flag,'iterations',iterations,'working_dimension',ell, ...
    'value_error',metrics.value_error,'residual',metrics.residual, ...
    'orthogonality',metrics.orthogonality,'pair_error',metrics.pair_error, ...
    'actual_backend',method,'used_builtin_svds',false,'fallback_reason','');
if isfield(metrics,'actual_backend'), row.actual_backend=metrics.actual_backend; row.used_builtin_svds=metrics.used_builtin_svds; row.fallback_reason=metrics.fallback_reason; end
end

function write_json(filename,data)
file = fopen(filename,'w');
if file < 0, error('benchmark:WriteFailed','Cannot write %s.',filename); end
cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
fprintf(file,'%s',jsonencode(data,'PrettyPrint',true));
end
