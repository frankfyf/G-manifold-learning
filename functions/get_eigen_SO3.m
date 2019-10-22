% function [ Evec, Eval ] = get_eigen_SO3( wigner_D, list, eigen_num, n )
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Build a series of affinity matrices based on the edge
% connections, then do eigen-decomposions. The transformation is SO(3)
% 
% Parameters : wigner_D          -- transformation on each edge 
%              list              -- edge connection list  
%              eigen_num         -- truncation of eigenvalues
%              n                 -- number of nodes
% 
% Return     : Eval              -- eigenvalues, a 1 by k_max cell
%              Evec              -- eigenvectors, a 1 by k_max cell

function [ Evec, Eval ] = get_eigen_SO3( wigner_D, list, eigen_num, n )

k_max = size(eigen_num,2); % maximum frequency

Evec = cell(1,k_max);
Eval = cell(1,k_max);
W = cell(1,k_max);
parfor k = 0:k_max-1
    W{1,k+1} = zeros((2*k+1)*n,(2*k+1)*n);
    for i = 1:size(list,1)
        W{1,k+1}((2*k+1)*(list(i,1)-1)+1: (2*k+1)*(list(i,1)-1)+2*k+1, ...
        (2*k+1)*(list(i,2)-1)+1: (2*k+1)*(list(i,2)-1)+2*k+1) = wigner_D{i,k+1};
    end
    tmp_D = sum(abs(sparse(list(:,1), list(:,2), 1, n, n)),2);
    D = repmat(tmp_D',2*k+1,1);
    D = D(:);
    W{1,k+1} = bsxfun(@times, 1./sqrt(D), W{1,k+1});
    W{1,k+1} = bsxfun(@times, 1./sqrt(D).', W{1,k+1});
    [ u, d ] = eigs(W{1,k+1}, eigen_num(k+1)*(2*k+1), 'lm');
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    sorted_eigvec = u(:, id);
    Evec{k+1} = sorted_eigvec;
    Eval{k+1} = sorted_eigval;
end

end

