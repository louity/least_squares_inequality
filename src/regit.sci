function [dist_v, x_v] = regit(A_m, b_v, alpha, N)
    // Search for the regularized least square inequalities solution.
    [m, n] = size(A_m)
    x_v = rand(n, 1);
    dist_v=[]
    
    for k = 1:N
        xold_v = x_v;
        x_v = -1/alpha * A_m' * max(A_m*x_v - b_v,0);
        dist_v = [dist_v sum(abs(x_v - xold_v))];
    end
endfunction

function [distend_v, index_v, pindex_v] = nbindex(A_m, b_v , alpha_v, N)
    [m,n] = size(A_m)
    distend_v = []
    index_v = []
    
    for alpha = alpha_v
        [dist_v, x_v] = regit(A_m, b_v, alpha, N)
        distend_v = [distend_v dist_v($)]
        index_v = [index_v sum(A_m*x_v - b_v <= 0)]
    end
    
    pindex_v = index_v / m
endfunction

function [] = plotindex(A_m, b_v, alpha_v, N)
    [distend_v, index_v, pindex_v] = nbindex(A_m, b_v , alpha_v, N)
    
    figure;
    clf;
    plot2d(alpha_v, distend_v)
    title(['Distance à la solution, N=', string(N)])
    a=gca();
    fig=a.children.children;
    fig.thickness=2;
    
    figure;
    clf;
    plot2d(alpha_v, pindex_v)
    title(['Pourcentage d indices corrects, N=', string(N)])
    a=gca();
    fig=a.children.children;
    fig.thickness=2;
endfunction

