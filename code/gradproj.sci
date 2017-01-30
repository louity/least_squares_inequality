function [xf_v, error_v, X_m] = gradproj(A_m, b_v, x0_v, N, rho)
    // Use projected gradient to solve the linear inequalities least-
    // square problem.
    // Input variables:
    //  A_m -- a matrix of size mxn.
    //  b_v -- a vector of size mx1.
    //  x0_v -- an input vector of size nx1.
    //  N -- number of iterations.
    //  rho -- step size of the algorithm.
    //
    // Output variables:
    //  xf_v -- a vector of size nx1
    //  error_v -- a vector of size Nx1
    //  X_m -- a matrix of size nxM storing all the values of x
    
    [m, n] = size(A_m);
    error_v = [];
    X_m = [x0_v];
    
    tildeA_m = [A_m -A_m eye(m,m)];
    
    Xf_v = [x0_v ; 0*x0_v ; 0*b_v];
    for k = 1:N
        Xint_v = Xf_v;
        Xf_v = Xf_v - 2*rho*tildeA_m'*(tildeA_m*Xf_v - b_v);
        Xf_v = (Xf_v>0) .* Xf_v;
        
        X_m = [X_m Xf_v(1:n)-Xf_v(n+1:2*n)];
        
        error_v = [error_v sum(abs(Xf_v-Xint_v))];
    end
    
    xf_v = Xf_v(1:n) - Xf_v(n+1:2*n)
endfunction

function [] = plot_gradproj(X_m, error_v, A_m, b_v)
    figure;
    clf;
    plot2d(1:length(error_v), log(abs(error_v)))
    a = gca();
    fig = a.children.children;
    fig.thickness = 2;
    ylabel("Erreur itérée logarithmique")
    xlabel("Itérations")
    
    res=[]
    for k = 1:size(X_m,2)
        res = [res A_m*(X_m(:,k))-b_v];
    end
    
    for k = 1:size(X_m,1)
        figure;
        clf;
        plot(1:size(X_m,2), res(k,:))
        a = gca();
        fig = a.children.children;
        fig.thickness = 2;
        ylabel("Valeur résiduelle")
        xlabel("Itérations")
    end
endfunction
