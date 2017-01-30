function [r] = isconsistent_lowrk(A_m, b_v)
    // Returns 1 if the system Ax<=b is consistent, -1 otherwise if 
    // the matrix is of low rank (0 or 1).
    // Input arguments:
    //      A_m -- a matrix of size n x m.
    //      b_v -- a vector of size n x 1
    //
    // Output arguments:
    //      r -- a real number equal to 1 if the system is consistent,
    //           -1 otherwise. r is equal to 0 if the matrix is of
    //           rank superior or equal to 2.
    
    [n, m] = size(A_m);
    [nb, mb] = size(b_v);

    if rank(A_m) >= 2 then
        disp('Error: rank(A)>=2');
        r = 0;
    elseif rank(A_m) == 0 then
        if sum(b_v >= 0) == nb then
            r = 1;
        else
            r = -1;
        end
    else
        k = 1;
        btest_v = A(:, k);
        nbtest = norm(btest_v);
        while (nbtest == 0) & (k <= m)
            btest_v = A(:, k+1);
            nbtest = norm(btest_v);
            k = k+1;
        end

        c = btest_v .* b_v;
        if or([sum(c >= 0) == n, sum(c <= 0) == n]) then
            r = 1;
        else 
            r = -1;
        end
    end

endfunction
