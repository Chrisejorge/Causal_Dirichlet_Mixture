function newLoglike = computeloglike_continous(T,X,Y,TX_matrix,X_matrix,p1,p2,z,k,kjs,kljs,beta_coeff,sigma2_coeff,gamma_coeff,x_pi_param,x_mean_param,x_var_param)
% Calculation of likelihood of data

loglikeregy = 0;
loglikeregt = 0;
loglikeregx = 0;
gamma_count = 1;
x_param_count = 1;
for jj = 1:k
    tmpidx1 = z(:,1) == jj;
    tx_tmp = TX_matrix(tmpidx1,:);
    y_tmp = Y(tmpidx1,:);
%     loglikeregy = loglikeregy + sum( sum( -log(sqrt(sigma2_coeff(jj,:))) -log(2*pi)/2 ...
%                     - 0.5 * ((y_tmp - tx_tmp * beta_coeff(jj,:)')./ sqrt(sigma2_coeff(jj,:))).^2 ) );       
    loglikeregy = loglikeregy + sum(sum(lognormpdf(y_tmp,tx_tmp * beta_coeff(jj,:)',sqrt(sigma2_coeff(jj,:))) ));

    kj = kjs(jj);
    for l = 1:kj
        gamma_coef_current = gamma_coeff(gamma_count,:);
        tmpidx2 = z(:,1) == jj & z(:,2) == l;
        x_tmp = X_matrix(tmpidx2,:);
        t_tmp = T(tmpidx2,:);
%         loglikeregt = loglikeregt + sum( log( binopdf(t_tmp,1,expit(x_tmp * gamma_coef_current')) + eps));
        loglikeregt = loglikeregt + sum( logbernpdf(t_tmp,expit(x_tmp * gamma_coef_current')) );
        gamma_count = gamma_count+1;
        klj = kljs(jj,l);
        
        logprodx1 = 0;
        logprodx2 = 0;
        for h = 1:klj
            tmpidx3 = z(:,1) == jj & z(:,2) == l & z(:,3) == h;
                        
            if p1 > 0
%                 logprodx1 = sum( sum( log(binopdf(X(tmpidx3,1:p1),ones(sum(tmpidx3),p1),repmat(x_pi_param(x_param_count,:),sum(tmpidx3),1) ) + eps) ) );
                 logprodx1 = sum(sum( logbernpdf(X(tmpidx3,1:p1),x_pi_param(x_param_count,:)) ));
            end
            
            if p2 > 0
%                 logprodx2 = sum( sum(-log(sqrt(x_var_param(x_param_count,:))) -log(2*pi)/2 ...
%                     - 0.5 * ((X(tmpidx3,p1+1:p1+p2)-x_mean_param(x_param_count,:))./sqrt(x_var_param(x_param_count,:))).^2) );
                  logprodx2 = sum(sum(lognormpdf(X(tmpidx3,p1+1:p1+p2),x_mean_param(x_param_count,:),sqrt(x_var_param(x_param_count,:))) ));
            end
            loglikeregx = loglikeregx + logprodx1 +  logprodx2;
            
            x_param_count = x_param_count + 1;
        end % for h
    end % for l
end % for jj
newLoglike = loglikeregy + loglikeregt + loglikeregx;
end
