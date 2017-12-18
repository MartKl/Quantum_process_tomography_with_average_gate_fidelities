% The goal is to demonstrate that MOSEK has a memory leak.

% In this example, a rank-1 matrix is recovered from dimension deficient data

% initialize
clear
% run /home/kliesch/.cvx/cvx_setup.m
% cvx_solver Mosek

% matrix dimension 
n = 15;
% number of tests
reps = 1;

for j=1:reps
    % draw rank-1 matrix and vectorize
    M0 = randn([n,1])*randn([n,1])';
    vecmat = reshape(M0,[n^2,1]);
    % draw linear map
    m = 6*n;
    A = randn(m,n^2);
    % y = incomplete information about M0
    y = A*vecmat ;%+ .01*randn([m,1]);
    %optimize
    cvx_begin sdp %quiet
        variable X(n,n)
        variable Y(n,n)
        variable M(n,n)
        % minimize trace norm of M
        minimize(trace(X+Y)/2)
        [(X+X')/2        M      ;
               M'        (Y+Y')/2] >= 0
        % subject to y == A*M
        y == A*reshape(M,[n^2,1])
    cvx_end
    
%     cvx_begin sdp quiet
%         variable M(n,n)
%         minimize(norm_nuc(M))
%         y == A*reshape(M,[n^2,1])
%     cvx_end

    
    
    %y = A*reshape(X0,[n^2,1])
   
    
    cvx_clear
    % compare M to M0
    fmat(j) = norm(reshape(M,[n^2,1]) - vecmat);
end

% largest deviation
largest = max(abs(fmat));
display(['largest deviation= ' num2str(largest)]);
