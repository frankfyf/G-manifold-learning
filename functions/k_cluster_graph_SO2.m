% function [ list, angle ] = k_cluster_graph_SO2( p, index )
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/17   
% 
% Description: Generate a random graph with the given random rewiring 
% probability and the corresponding SO(2) transformation between connected
% nodes
% 
% Parameters : p                 -- probability of random rewiring model
%              index             -- K length cell, each cell contains a
%                                   vector recording the index of data
%                                   points in each cluster
% 
% Return     : list              -- edge connection list
%              angle             -- transformation (rotation angle) on each
%                                   edge

function [ list, angle ] = k_cluster_graph_SO2(p, index)

% Get parameters 
K = numel(index); % number of clusters
n_points = zeros(1,K);
for i = 1:K
    n_points(i) = numel(index{i}); 
end
n = sum(n_points); 

% Generate random transformation in SO(2)
r = 2*pi*rand(1,n);

% build the list of edge connection
list = zeros(sum(n_points.*(n_points-1)/2), 2);
angle = zeros(sum(n_points.*(n_points-1)/2), 1);
count = 1;
for kk = 1:K
    for i = 1:n_points(kk)
        for j = i+1:n_points(kk)
            if(rand(1) <= p) % with probability p we keep the edge connection
                list(count,:) = [index{1,kk}(i), index{1,kk}(j)];
                angle(count,1) = r(index{1,kk}(i)) - r(index{1,kk}(j)); 
            else
                tmp_num = randi(n);
                if(rand(1) < 0.5)
                    while(tmp_num == index{1,kk}(i))
                        tmp_num = randi(n);
                    end
                    list(count,:) = [tmp_num,index{1,kk}(j)];
                    angle(count,1) = 2*pi*rand(1);
                else
                    while(tmp_num == index{1,kk}(j))
                        tmp_num = randi(n);
                    end
                    list(count,:) = [index{1,kk}(i), tmp_num];
                    angle(count,1) = 2*pi*rand(1);
                end
            end
        count = count+1;
        end
    end
end

% Remove the same edge connection
[list, id] = sort(list,2, 'ascend');
for i = 1:size(id,1)
    if id(i,1) == 2
        angle(i,1) = -angle(i,1);
    end
end
[list, id, ~] = unique(list,'rows');
angle = angle(id,:);
list = [list; list(:,end:-1:1)];
angle = [angle; -angle];

% Add self connection
tmp_list = repmat([1:n]',1,2);
tmp_angle = zeros(n,1);
list = [list;tmp_list];
angle = [angle; tmp_angle];

end

