% function [ affinity ] = aff_opt(Evec, Eval, t)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Compute the optimal alignment affinity between nodes, based on 
% the eigen-decomposition. The transformation is SO(2). The alignment is
% computed by 
% 
% Parameters : Evec              -- eigenvectors, a 1 by k_max cell
%              Eval              -- eigencalues, a 1 by k_max cell
%              t                 -- diffusion time
% 
% Return     : affinity          -- a n by n matrix, the affinity between n nodes.

function [ affinity ] = aff_opt(Evec, Eval, t)

% Build a n by n by k_max tensor, the alignment between (i,j) pair is
% computed by the FFT of the (i,j,:)th slice 
k_max = size(Evec,2);
n = size(Evec{1}, 1);
align_tensor = zeros(n,n, k_max);

% Weighted matrices filtering.
for i = 1:k_max
    evec = Evec{i};
    evec_norm = sqrt(sum(abs(evec).^2, 2));
    evec = bsxfun(@times, evec, 1./evec_norm);
    align_tensor(:,:,i) = evec*diag(Eval{i}.^(2*t))*(evec');
end

% FFT and compute the optimal alignment affinity
affinity = zeros(n,n);
alignment = zeros(n,n);
N = 1024; % number of FFT points
tmp_affinity = zeros(1,N);
for i = 1:n
    for j = i+1:n
        tmp_affinity(1,1:k_max) = align_tensor(i,j,:);
        tmp_fft = abs(fft(tmp_affinity, N));
        affinity(i,j) = max(tmp_fft);
        affinity(j,i) = affinity(i,j);
        alignment(i,j) = 360*(find(tmp_fft == affinity(i,j),1)-1)/N;
        if(alignment(i,j) > 180)
            alignment(i,j) = alignment(i,j) - 360;
        end
        alignment(j,i) = -alignment(i,j);
    end
end

end

