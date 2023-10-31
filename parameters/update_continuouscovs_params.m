function [newmean,newvar] = update_continuouscovs_params(X,p1,p2,z,k,kjs,kljs,nu0,tau0,c0,mu0)
% Update continuous covariate parameters
% Within each cluster parameters are stored in a long vector nXClus x 1
% Continuous covariates have 2 parameters (newmean, newvar)

% Initialize vectors nXClus is the total number of clusters (t,x and y clusters)
nXClus = size(unique(z,'rows'),1);

if p2 > 0
    newmean = zeros(nXClus,  p2);
    newvar = zeros(nXClus, p2);
    
    % beta(a0,b0) prior,a0 = 1, b0 = 1
    count = 1;
    for jj = 1:k
        kj = kjs(jj);
        for l = 1:kj
            klj = kljs(jj,l);
            for h = 1:klj
                tmpidx = z(:,1) == jj & z(:,2) == l & z(:,3) == h;
                
                x_tmp = X(tmpidx,p1+1:p1+p2);
                % posterior is beta
                newvar(count,:) = update_covs_variance(nu0,tau0,c0,mu0,x_tmp);
                
                newmean(count,:) = update_covs_mean(c0,mu0,x_tmp,newvar(count,:));
                
                count = count + 1;
            end % for h
        end % for l
    end % for jj
    
else
    error('Requires at least one dimension (p2 > 0).');
end

% end