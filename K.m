function out = K(m,n)
% Commutation matrix: vec(A') = K(m,n)*vec(A) for an mxn matrix A
    I = reshape(1:m*n, [m, n]); % initialize a matrix of indices of size(A)
    I = I'; % Transpose it
    I = I(:); % vectorize the required indices
    out = eye(m*n); % Initialize an identity matrix
    out = out(I,:); % Re-arrange the rows of the identity matrix