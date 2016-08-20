function X = pinv_sparse(A,tol)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS(class(A)).
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%
%   Class support for input A: 
%      float: double, single
%
%   See also RANK.
 
%   Copyright 1984-2013 The MathWorks, Inc. 

[U,S,V] = svds(A,size(A,2));
s = diag(S);
if nargin < 2 
    tol = max(size(A)) * eps(norm(s,inf));
end
r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./s(:);
X = sparse(bsxfun(@times,V,s.'))*U';
