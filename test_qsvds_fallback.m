function test_qsvds_fallback
%TEST_QSVDS_FALLBACK Force native errors/poor factors; verify pure recovery.
rng(919,'twister');
A=quaternion_to_complex(randn(32,25)+1i*randn(32,25),randn(32,25)+1i*randn(32,25));
reference=svd(A); expected=reference(1:2:4);
folder=tempname; mkdir(folder); addpath(folder,'-begin');
cleanup=onCleanup(@()cleanup_guard(folder)); %#ok<NASGU>
saved_warning=warning('off','qsvds:NoConvergence');
warning_cleanup=onCleanup(@()warning(saved_warning)); %#ok<NASGU>
for mode=1:2
    file=fopen(fullfile(folder,'svds.m'),'w');
    if mode==1
        fprintf(file,'function varargout=svds(varargin)\nerror(''qsvds:InjectedSVDSFailure'',''Test injection.'');\nend\n');
    else
        fprintf(file,['function [U,S,V,flag]=svds(A,k,varargin)\n' ...
            '[V,~]=qr(randn(size(A,2),k)+1i*randn(size(A,2),k),0);\n' ...
            'U=zeros(size(A,1),k); S=eye(k); flag=0;\nend\n']);
    end
    fclose(file); rehash; clear qsvds_native svds
    assert(strcmp(which('svds'),fullfile(folder,'svds.m')));
    state=rng;
    [U,E,V,info]=QSVDS(A,32,25,2);
    assert(isequal(state,rng));
    assert(info.flag==0 && info.used_builtin_svds && ~isempty(info.fallback_reason));
    assert(~isempty(strfind(info.method,'krylov'))); %#ok<STREMP>
    assert(norm(diag(E(1:2,1:2))-expected)<1e-8*expected(1));
    assert(norm(A*V-U*E,'fro')/norm(A,'fro')<1e-9);
end
fprintf('Native-exception and bad-native-factor fallbacks passed; no full decomposition requested.\n');
end
function cleanup_guard(folder)
rmpath(folder); clear qsvds_native svds
delete(fullfile(folder,'svds.m')); rmdir(folder);
end
