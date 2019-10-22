% function [ affinity ] = aff_ps_SO3( Evec, Eval, t )
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/19   
% 
% Description: Compute the power spectrum affinity between nodes, based on 
% the eigen-decomposition. The transformation is SO(3). 
% 
% Parameters : Evec              -- eigenvectors, a 1 by k_max+1 cell
%              Eval              -- eigencalues, a 1 by k_max+1 cell
%              t                 -- diffusion time
% 
% Return     : affinity          -- a n by n matrix, the affinity between n nodes.

function [ affinity ] = aff_ps_SO3( Evec, Eval, t )

% Weight matrices filtering
k_max = size(Evec, 2);
n = size(Evec{1}, 1);
W = cell(1,k_max);
for i = 0:k_max-1
    eval = Eval{i+1}.^(2*t);
    evec = Evec{i+1};
    evec_norm = sqrt(sum(abs(evec).^2, 2));
    evec = bsxfun(@times, evec, 1./evec_norm);
    W{i+1} = bsxfun(@times, evec, eval.')*evec'; 
end

% Compute the power spectrum affinity
affinity = zeros(n,n,k_max);
for i = 1:n
    for j = i+1:n
        for kk = 1:k_max-1
            tmp = W{kk+1}((i-1)*(2*kk+1)+1:i*(2*kk+1), (j-1)*(2*kk+1)+1:j*(2*kk+1));
            affinity(i,j,kk+1) = trace(tmp*tmp');
        end
        disp(['Finish the power spectrum affinity for nodes ', num2str(i), ' and ', num2str(j)]);
    end
end

affinity = sum(affinity,3)/(k_max-1);
affinity = affinity + affinity';

end

