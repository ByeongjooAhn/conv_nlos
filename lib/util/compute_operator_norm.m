function [L1, L2]= compute_operator_norm(A, AS, sx)

    % computes the operator norm for a linear operator AS on images with size sx, 
    % which is the square root of the largest eigenvector of AS*A.
    % http://mathworld.wolfram.com/OperatorNorm.html
    vec_size = prod(sx);

    %Compute largest eigenvalue (in this case arnoldi, since matlab
    %implementation faster than power iteration)
    opts.tol = 1.e-3;
    opts.maxit = 50;
    lambda_largest = eigs(@(x)ASAfun(x, A, AS, sx), vec_size, 1,'LM', opts);
    L1 = sqrt(lambda_largest);
    lambda_smallest = eigs(@(x)ASAfun(x, A, AS, sx), vec_size, 1,'SM', opts);
    L2 = sqrt(lambda_smallest);
return;

function ASAx = ASAfun(x, A, AS, sx)
    x_img = reshape(x,sx);
    ASAx = AS(A(x_img));
    ASAx = ASAx(:);
return;