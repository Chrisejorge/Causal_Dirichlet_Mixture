function newpi = update_binarycovs_params(X,p1,z,k,kjs,kljs,a0,b0)
% Update binary covariate parameters
% Within each cluster parameters are stored in a long vector nXClus x 1
% Binary covariates have 1 parameter (x_pi_param)
% Initialize vectors nXClus is the total number of clusters (t,x and y clusters)
nXClus = size(unique(z,'rows'),1);
if p1 > 0
    newpi = zeros(nXClus, p1);
    
    %  Update parameter for binary covariates
    
    % beta(a0,b0) prior, a0 = 1, b0 = 1
    count = 1;
    for jj = 1:k
        kj = kjs(jj);
        for l = 1:kj
            klj = kljs(jj,l);
            for h = 1:klj
                tmpidx = z(:,1) == jj & z(:,2) == l & z(:,3) == h;
                x_tmp = X(tmpidx,1:p1);
                % posterior is beta
                newpi(count,:) = betarnd(a0 + sum(x_tmp,1),b0 + size(x_tmp,1)- sum(x_tmp,1));
                count = count + 1;
            end % for h
        end % for l
    end % for jj
    
else
    error('Requires at least one dimension (p1 > 0).');
end

% end