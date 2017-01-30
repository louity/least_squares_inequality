function [dist_v, xf_v] = han(A_m, b_v, N)
    [m,n] = size(A_m);
    xf_v = rand(n,1);
    dist_v = [];
    
    for k = 1:N
        AI_m = A_m(A_m*xf_v-b_v>=0, :);
        bI_v = b_v(A_m*xf_v-b_v>=0);
        
        xold_v = xf_v
        xf_v = xf_v - pinv(AI_m)*(AI_m*xf_v-bI_v)
        dist_v = [dist_v sum(abs(xf_v-xold_v))]
    end
endfunction

function [nbit] = nbit_han(n_v, nit, Nmax)
    nbit= [];
    
    for n=n_v
        nbitint = [];
        
        for j=1:nit
        A_m = rand(n,n);
        b_v = rand(n,1);
        [dist_v, xf_v] = han(A_m, b_v, Nmax);
        nbitint = [nbitint sum(1.*(dist_v>0))];
        end
        
        nbit = [nbit mean(nbitint)];
    end
endfunction