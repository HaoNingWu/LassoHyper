function beta = l2_beta(w,A,f,lambda,mu)
% For coefficients of Lasso hyperinterpolation
    Af = A'*diag(w)*f;
    beta = Af./(1+lambda*mu.^2);
end
