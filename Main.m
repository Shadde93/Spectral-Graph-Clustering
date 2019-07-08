clc 
clear all
%% Read file and choose k

Example1 = csvread('example1.dat');
Example2 = csvread('example2.dat');

%% Get affinity matrix and plot the graphical model, choose which example

[G, A] = get_A_G(Example1);

figure(1)
p = plot(G,'-ok','Layout','force');
title('Cluster of nodes');

%% Spectral clustering algorithm with normalized laplacian matrix, choose k here

k = 2;

[L,idx] = S_C_Alg_Lnorm(A,k);

% set cluster points to the nodes on the prevoius graph

color = hsv(k);

for i=1:k
    cluster_points = find(idx==i);
    highlight(p,cluster_points,'NodeColor',color(i,:))
end
title('Spectral clustered nodes');

% Plot affininty patterns and fiddler vectors
figure(2)
spy(A)
title('Sparsity patterns');

figure(3)
[L_v,~] = eigs(L, k, 'LM'); %change SA to LM

fiedler_vector = sort(L_v(:,k));
plot(fiedler_vector,'-o');
title('Sorted Fiedler Vector');
%% Spectral clustering algorithm with laplacian matrix, choose k here

k = 4;

[L,idx] = S_C_Alg_L(A,k);

% set cluster points to the nodes on the prevoius graph

color = hsv(k);

for i=1:k
    cluster_points = find(idx==i);
    highlight(p,cluster_points,'NodeColor',color(i,:))
end
title('Spectral clustered nodes');

% Plot affininty patterns and fiddler vectors

figure(2)
spy(A)
title('Sparsity patterns');

figure(3)
[L_v,~] = eigs(L, k, 'SA');

fiedler_vector = sort(L_v(:,k));
plot(fiedler_vector,'-o');
title('Sorted Fiedler Vector');

%% Functions

function [graph_model, affinity_matrix] = get_A_G(edges)

max_node = max(max(edges)); 
u = edges(:,1);
v = edges(:,2);
A_sparse = sparse(u, v, 1, max_node, max_node); 

af = full(A_sparse);
affinity_matrix = full(spones(af));
 %ignores the diagonal entries of the adjacency matrix A and does not add self-loops to the graph.
graph_model = graph(affinity_matrix,'OmitSelfLoops');

end


function [Laplacian_matrix, cluster_list] = S_C_Alg_Lnorm(A_matrix,k)

D = diag(sum(A_matrix,2));
%D diagonal mateix with the eigenvalues, v is the eigen vectors (independant columns)
Laplacian_matrix = (D^(-1/2))*A_matrix*(D^(-1/2)); 
%largest magnitude LM
[X,~] = eigs(Laplacian_matrix, k, 'LM'); 

Y = X./sqrt(sum(X.^2,2));

cluster_list = kmeans(Y,k);
end

function [Laplacian_matrix, cluster_list] = S_C_Alg_L(A_matrix,k)

D = diag(sum(A_matrix,2));
%D diagonal mateix with the eigenvalues, v is the eigen vectors (independant columns)
Laplacian_matrix = D-A_matrix; 
%largest magnitude LM
[X,~] = eigs(Laplacian_matrix, k, 'SA'); 

Y = X./sqrt(sum(X.^2,2));

cluster_list = kmeans(Y,k);
end