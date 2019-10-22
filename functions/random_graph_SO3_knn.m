% function [ list, wigner_D ] = random_graph_SO3_knn(dist, kappa, p, k_max)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/21   
% 
% Description: Generate a random graph G = (V,E) by connecting kappa-nearest 
% neighbor and apply the random rewiring model for adding noise. The 
% SO(3) transformation between nodes is included. 
% 
% Parameters : dist              -- n by n matrices, distance between n nodes  
%              kappa             -- kappa-nearest neighbor
%              p                 -- probability of random rewiring model  
%              k_max                 -- maximum frequency
% 
% Return     : list              -- A |E| by 2 list of edges
%              wigner_D          -- The SO(3) transformation on each edge

function [ list, wigner_D ] = random_graph_SO3_knn(dist, kappa, p, k_max)

n = size(dist,1); % number of data points

% Assign SO(3) transformation for each node
angle = quat2eul(qrand(n)', 'ZYZ');
angle = mod(angle,2*pi);
wigner_l = cell(n,k_max+1);
for i = 1:n
    for k = 0:k_max
        [wigner_l{i,k+1},~,~] = wignerD(k,angle(i,1), angle(i,2), angle(i,3));
    end
end

% Find the kappa nearest neighbors list
dist = dist - diag(diag(dist));
[~,id] = sort(dist,2, 'descend');
id = id(:,1:kappa);
list = [repmat([1:n]', kappa,1), id(:)];
list = unique(sort(list,2,'ascend'), 'rows');

% Compute the alignment and do rewiring for each edge
len = size(list,1);
wigner_D = cell(2*len, k_max+1);
for i = 1:len
    if rand(1) <= p % With probability p, we keep the connection
        for k = 0:k_max
            wigner_D{i,k+1} = wigner_l{list(i,1),k+1}*wigner_l{list(i,2),k+1}';
            wigner_D{i+len,k+1} = wigner_D{i+len,k+1}';
        end
    else % With probability 1-p, we random rewiring the edge and transformation
        tmp_angle = quat2eul(qrand(1)', 'ZYZ');
        tmp_angle = mod(tmp_angle,2*pi);
        for k = 0:k_max
            [wigner_D{i,k+1},~,~] = wignerD(k,tmp_angle(1,1), tmp_angle(1,2), tmp_angle(1,3));
            wigner_D{i+len,k+1} = wigner_D{i+len,k+1}';
        end
        tmp_id = randi(n);
        if rand(1) < 0.5
            while (tmp_id == list(i,1) || tmp_id == list(i,2))
                tmp_id = randi(n);
            end
            list(i,1) = tmp_id;
        else
            while (tmp_id == list(i,1) || tmp_id == list(i,2))
                tmp_id = randi(n);
            end
            list(i,2) = tmp_id;
        end
    end
end

list = [list; list(:,[2,1])];

% Remove the same connection
[list, id] = sort(list,2, 'ascend');
for i = 1:size(id,1)
    if id(i,1) == 2
        for k = 0:k_max
            wigner_D{i,k+1} = wigner_D{i,k+1}';
        end
    end
end
[list, id, ~] = unique(list,'rows');
wigner_D = wigner_D(id,:);

len = size(list,1);
list = [list; list(:,end:-1:1)];
tmp_wigner = cell(len,k_max+1);
for i = 1:len
    for k = 0:k_max
        tmp_wigner{i,k+1} = wigner_D{i,k+1}';
    end
end
wigner_D = [wigner_D; tmp_wigner];

end

