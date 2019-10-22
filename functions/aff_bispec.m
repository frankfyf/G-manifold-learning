% function [ affinity ] = aff_bispec(Evec, Eval, t)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/19   
% 
% Description: Compute the Bispectrum affinity between nodes, based on 
% the eigen-decomposition. The transformation is SO(2)
% 
% Parameters : Evec              -- eigenvectors, a 1 by k_max cell
%              Eval              -- eigencalues, a 1 by k_max cell
%              t                 -- diffusion time
% 
% Return     : affinity          -- a n by n matrix, the affinity between n nodes.

function [ affinity ] = aff_bispec(Evec, Eval, t)

k_max = size(Evec, 2);
n = size(Evec{1}, 1);
W = zeros(n, n, k_max + 1);

% Weighted matrices filtering 
for i = 1:k_max
    eval = Eval{i}.^(2*t);
    evec = Evec{i};
    evec_norm = sqrt(sum(abs(evec).^2, 2));
    evec = bsxfun(@times, evec, 1./evec_norm);
    W(:, :, i) = bsxfun(@times, evec, eval.')*evec'; 
end

% Compute the bispectrum affinity
affinity = zeros(n, n);
for k2 = 1:k_max
    for k1 = 1:k2 - 1
        affinity = affinity + W(:, :, k1).*conj(W(:, :, k2)).*W(:, :, k2 - k1);
    end
end

end

