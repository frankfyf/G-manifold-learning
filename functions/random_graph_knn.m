% function [ list, angle ] = random_graph_knn(dist, kappa, p, q)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/16   
% 
% Description: Generate a random graph G = (V,E) by connecting kappa-nearest 
% neighbor and apply the random rewiring model for adding noise. The 
% SO(2) transformation between nodes is included. 
% 
% Parameters : dist              -- n by n matrices, distance between n nodes  
%              kappa             -- kappa-nearest neighbor
%              p                 -- probability of random rewiring model  
%              q                 -- quaternions
% 
% Return     : list              -- A |E| by 2 list of edges
%              angle             -- A |E| by 1 list of corresponding transformations

function [ list, angle ] = random_graph_knn(dist, kappa, p, q)

n = size(dist,1); % number of data points

% Find the kappa-nearest neighbors list 
dist = dist - diag(diag(dist));
[~,id] = sort(dist,2, 'descend');
id = id(:,1:kappa);
list = [repmat([1:n]', kappa,1), id(:)];
list = unique(sort(list,2,'ascend'), 'rows'); % remove repeated rows

% Compute the transformation and apply random rewiring model
angle = zeros(size(list,1),1);
for i = 1:size(list,1)
    if rand(1) <= p % With probability p, we keep the right edge connection and transformations
        angle(i) = pi*q_to_inplanerot(list(i,1), list(i,2), q)/180;
    else % With probability 1-p, random rewire the edge connection and transformation
        angle(i) = 2*pi*rand(1);
        tmp_id = randi(n);
        if rand(1) < 0.5
            while tmp_id == list(i,1)
                tmp_id = randi(n);
            end
            list(i,1) = tmp_id;
        else
            while tmp_id == list(i,2)
                tmp_id = randi(n);
            end
            list(i,2) = tmp_id;
        end
    end
end
list = [list; list(:,[2,1])];
angle = [angle,-angle];
            
end

