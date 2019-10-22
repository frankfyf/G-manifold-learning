% function [ cg_m ] = cg_matrix( j1, j2, k_max )
% 
% Author:   Yifeng Fan (yifengf2@illinois.edu)
% Date:     2019/10/21   
% 
% Description: Compute the Clebsch-Gordan matrix with frequency j1, j2 and
% j. Each entry is a Clebsch-Gordan coefficient, and the output matrix is 
% (2*j1+1)*(2*j2+1) by n, where n = \sum_{i = |j1-j2|}^{min{j1+j2, k_max}}.
% 
% Parameters : j1                -- frequency j1
%              j2                -- frequency j2
%              k_max             -- maximum frequency
% 
% Return     : cg_m              -- Clebsch-Gordan matrix

function [ cg_m ] = cg_matrix( j1, j2, k_max )

cg_m = cell(1, max(min(j1+j2,k_max)+1 - abs(j2-j1), 0));
for j = abs(j2-j1):min(j1+j2,k_max)
    tmp = zeros((2*j1+1)*(2*j2+1), 2*j+1);
    for m1 = -j1:j1
        for m2 = -j2:j2
            for m = -j:j
                if( m == m1 + m2 )
                    tmp((m1+j1)*(2*j2+1)+m2+j2+1, m+j+1) = ClebschGordan(j1,j2,j,m1,m2,m);
                end
            end
        end
    end
    cg_m{1,j - abs(j2-j1) + 1} = tmp;
end
cg_m = cell2mat(cg_m);

end

