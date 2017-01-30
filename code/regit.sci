function [dist_v, x_v] = fixedpoint_reg(A_m, b_v, alpha, N)
    // Search for the regularized least square inequalities solution.
    [m, n] = size(A_m)
    x_v = rand(n, 1);
    dist_v = []

    for k = 1:N
        xold_v = x_v;
        x_v = -1/alpha * A_m' * max(A_m*x_v - b_v,0);
        dist_v = [dist_v sum(abs(x_v - xold_v))];
    end
endfunction

function [dist_v, x_v] = gradient_descent_reg(A_m, b_v, alpha, N, step)
    //Gradient descent with fixed stepsize.

    [m, n] = size(A_m)
    x_v = rand(n, 1);
    dist_v=[]

    for k = 1:N
        xold_v = x_v;
        x_v = x_v - step*(2*A_m'*max(A_m*x_v - b_v,0) + 2*alpha*x_v);
        dist_v = [dist_v sum(abs(x_v - xold_v))];
    end
endfunction

function [distend_v, index_v, pindex_v] = nbindex(A_m, b_v , alpha_v, N, method, step)
    [m,n] = size(A_m)
    distend_v = []
    index_v = []

    for alpha = alpha_v
        if method == 1 then 
            [dist_v, x_v] = fixedpoint_reg(A_m, b_v, alpha, N)
        elseif method == 2 then
            [dist_v, x_v] = gradient_descent_reg(A_m, b_v, alpha, N, step)
        elseif method == 3 then
            [dist_v, x_v] = han(A_m, b_v, N)
        end

        distend_v = [distend_v dist_v($)]
        index_v = [index_v sum(A_m*x_v - b_v <= 0)]
    end

    pindex_v = index_v / m
endfunction

function [] = plotindex_reg(A_m, b_v, alpha_v, N, step)
    [distendfp_v, indexfp_v, pindexfp_v] = nbindex_reg(A_m, b_v , alpha_v, N, 1, step)
    [distendgrad_v, indexgrad_v, pindexgrad_v] = nbindex_reg(A_m, b_v , alpha_v, N, 2, step)

    figure;
    clf;
    plot2d(alpha_v, distendfp_v)
    plot2d(alpha_v, distendgrad_v)
    title(['Distance à la solution, N=', string(N)])
    legend(['fixed point' ; 'gradient descent'])
    a=gca();
    fig2=a.children(3).children;
    fig2.thickness=2;
    fig1=a.children(2).children;
    fig1.thickness=2;
    fig1.foreground=5;

    figure;
    clf;
    plot2d(alpha_v, pindexfp_v)
    plot2d(alpha_v, pindexgrad_v)
    title(['Pourcentage d indices corrects, N=', string(N)])
    legend(['fixed points' ; 'gradient descent'])
    a=gca();
    fig2=a.children(3).children;
    fig2.thickness=2;
    fig1=a.children(2).children;
    fig1.thickness=2;
    fig1.foreground=5
endfunction



