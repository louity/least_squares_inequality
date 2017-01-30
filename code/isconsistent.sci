function [r] = isconsistent(A_m, b_v)
    // Returns 1 if the system Ax<=b is consistent, -1 otherwise.
    // Input arguments:
    //      A_m -- a matrix of size n x m.
    //      b_v -- a vector of size n x 1
    //
    // Output arguments:
    //      r -- a real number equal to 1 if the system is consistent,
    //           0 otherwise.
    
    r = 0
    N = 0
    if rank(A_m)<2 then
        r = isconsistent_lowrk(A_m, b_v);
    end

    while r == 0
        disp(N)
        disp(size(A_m))
        [n, m] = size(A_m);
        Ait_m = -[A_m b_v]';
        bit_v = zeros(m+1,1);
        bit_v($)=1;

        [xsol_v, kerAit_m] = linsolve(Ait_m, -bit_v);

        if isempty(xsol_v)
            r = (-1)^N;
        else 
            A_m = -kerAit_m;
            b_v = xsol_v;
        end

        if rank(A_m)<2 then 
            r = isconsistent_lowrk(A_m, b_v);
            r = (-1)^(N-1) * r;
        end

        N=N+1;
    end
    
endfunction
