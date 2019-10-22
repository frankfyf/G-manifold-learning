% function [ list, wigner_D ] = k_cluster_graph_SO3( p, index, k_max)
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/19   
% 
% Description: Generate a random graph with the given random rewiring 
% probability and the corresponding SO(3) transformation between connected
% nodes. Note that here the transformation is expressed by wigner-D
% matrices.
% 
% Parameters : p                 -- probability of random rewiring model
%              index             -- K length cell, each cell contains a
%                                   vector recording the index of data
%                                   points in each cluster
%              k_max             -- maximum frequency
%               
% 
% Return     : list              -- edge connection list
%              wigner_D          -- The SO(3) transformation on each edge

function [ list, wigner_D ] = k_cluster_graph_SO3( p, index, k_max )

% Get parameters 
K = numel(index); % number of clusters
n_points = zeros(1,K);
for i = 1:K
    n_points(i) = numel(index{i}); 
end
n = sum(n_points); % the total number of nodes

% Assign the wigner_d matrix for each node
angle = quat2eul(qrand(n)', 'ZYZ'); 
angle = mod(angle,2*pi);
wigner_l = cell(n,k_max+1);
for i = 1:n
    for k = 0:k_max
        [wigner_l{i,k+1},~,~] = wignerD(k,angle(i,1), angle(i,2), angle(i,3));
    end
end

% Build the list of connection
list = zeros(sum(n_points.*(n_points-1)/2), 2);
wigner_D = cell(2*sum(n_points.*(n_points-1)/2), k_max+1);
count = 1;
for kk = 1:K
    for i = 1:n_points(kk)
        for j = i+1:n_points(kk)
            if(rand(1) <= p) % With probability p, we keep the connection
                list(count,:) = [index{1,kk}(i), index{1,kk}(j)];
                for k = 0:k_max
                    wigner_D{count,k+1} = wigner_l{index{1,kk}(i),k+1}*wigner_l{index{1,kk}(j),k+1}';
                    wigner_D{count+sum(n_points.*(n_points-1)/2),k+1} = wigner_D{count,k+1}';
                end
            else % With probability 1-p, we random rewiring the edge and transformation
                tmp_num = randi(n);
                if(rand(1) < 0.5)
                    while(tmp_num == index{1,kk}(i) || tmp_num == index{1,kk}(j))
                        tmp_num = randi(n);
                    end
                    list(count,:) = [tmp_num,index{1,kk}(j)];
                else
                    while(tmp_num == index{1,kk}(i) || tmp_num == index{1,kk}(j))
                        tmp_num = randi(n);
                    end
                    list(count,:) = [index{1,kk}(i), tmp_num];
                end
                tmp_angle = quat2eul(qrand(1)', 'ZYZ');
                tmp_angle = mod(tmp_angle,2*pi);
                for k = 0:k_max
                    [wigner_D{count,k+1},~,~] = wignerD(k,tmp_angle(1,1), tmp_angle(1,2), tmp_angle(1,3));
                    wigner_D{count+sum(n_points.*(n_points-1)/2),k+1} = wigner_D{count,k+1}';
                end
            end
        count = count+1;
        end
    end
end
list = [list; list(:,end:-1:1)];

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

% Add the self connection in SO(3)
tmp_list = zeros(n,2);
tmp_wigner = cell(n,k_max+1);
count = 1;
for kk = 1:K
    for i = 1:n_points(kk)
        tmp_list(count,:) = [index{1,kk}(i), index{1,kk}(i)];
        for k = 0:k_max
            tmp_wigner{count,k+1} = eye(2*k+1);
        end
        count = count+1;
    end
end
list = [list;tmp_list];
wigner_D = [wigner_D; tmp_wigner];

end

