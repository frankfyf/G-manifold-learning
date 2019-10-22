% function [ affinity ] = aff_ps(Evec, Eval, t)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Compute the power spectrum affinity between nodes, based on 
% the eigen-decomposition. The transformation is SO(2). Code is built based on
% the Multi-frequency vector diffusion maps
% (http://proceedings.mlr.press/v97/fan19a/fan19a.pdf)
% 
% Parameters : Evec              -- eigenvectors, a 1 by k_max cell
%              Eval              -- eigencalues, a 1 by k_max cell
%              t                 -- diffusion time
% 
% Return     : affinity          -- a n by n matrix, the affinity between n nodes.

function [ affinity ] = aff_ps(Evec, Eval, t)

k_max = size(Evec, 2);
n = size(Evec{1}, 1);

% Weight matrices filtering
W = zeros(n, n, k_max);
for i = 1:k_max
    eval = Eval{i}.^(2*t);
    %heval = 2*eval - eval.^2;
    heval = eval;
    evec = Evec{i};
    evec_norm = sqrt(sum(abs(evec).^2, 2));
    evec = bsxfun(@times, evec, 1./evec_norm);
    W(:, :, i) = bsxfun(@times, evec, heval.')*evec'; %denoising the matrices
end

affinity = zeros(n, n);
for k1 = 1:k_max
    affinity = affinity + abs(W(:, :, k1)).^2;
end

end

