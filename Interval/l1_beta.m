function beta = l1_beta(w,A,f,lambda,mu)
% For coefficients of Lasso hyperinterpolation
% 
    Af = A'*diag(w)*f;
%     for i = 0:L
%         if 2*Af(i+1)>lambda*mu(i+1)
%             beta(i+1) = (2*Af(i+1)-lambda*mu(i+1))/(2*1);
%         elseif 2*Af(i+1)<-lambda*mu(i+1)
%             beta(i+1) = (2*Af(i+1)+lambda*mu(i+1))/(2*1);
%         else
%             beta(i+1) = 0;
%         end
%     end
      beta = shrink(Af,lambda*mu);


end


function z = shrink(x, r)
    z = sign(x).*max(abs(x)-r,0);
end
