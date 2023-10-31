function [X,TX,X_matrix,TX_matrix] = update_missing(T,X,Y,TX_matrix,pt,p1,p2,X_mar,z,k,kjs,kljs,beta_coeff,x_pi_param,x_mean_param,x_var_param)
% Update missing coviriates
% unique clusters memberships
    
    xClusNum = NaN(size(X,1),1);
    count = 1;
    for jj = 1:k
        kj = kjs(jj);
        for l = 1:kj
            klj = kljs(jj,l);
            for h = 1:klj
                tmpidx = z(:,1) == jj & z(:,2) == l & z(:,3) == h;
                xClusNum(tmpidx) = count;
                count = count + 1;
            end % for h
        end % for l
    end % for jj
    
    % Binary covariates: assume binomial distribution
    if p1 > 0
        for r = 1:p1
            % Find missing observations
            tmpidx = X_mar(:,r) == 1;
            tmpx1 = TX_matrix(tmpidx,:);
            
            % Calculate probability that missing binary covariate is 1.
            % Set missing covariate to 0 and 1.
            tmpx0 = tmpx1;
            tmpx1(:, r + 1 + pt) = 1;
            tmpx0(:, r + 1 + pt) = 0;
            tmpexp1 = expit(sum(tmpx1 .* beta_coeff(z(tmpidx,1),:),2));
            tmpexp0 = expit(sum(tmpx0 .* beta_coeff(z(tmpidx,1),:),2));
            piidx = xClusNum(tmpidx);          
            tmppis = x_pi_param(piidx,r);
            tmpvalue1 = tmppis .* tmpexp1;
            tmpvalue0 = (1 - tmppis) .* tmpexp0;            
            tmprob =  tmpvalue1 ./ (tmpvalue1 + tmpvalue0);
            
            % Draw missing values from binomial distribution
            X(tmpidx,r) =  binornd(1, tmprob ,size(tmprob,1),1);
        end % for r
    end
    
   % Continuous covariates are updated using Metropolis Hastings
    if p2 > 0
        for r = p1+1:p1+p2
            tmpidx = X_mar(:,r) == 1;
            proposed = normrnd(X(tmpidx,r),0.1);
            propTX_matrix = TX_matrix(tmpidx,:);
            propTX_matrix(:,1+pt+r) = proposed;
            
            varidx = xClusNum(tmpidx);          
            loglikeprop = -log( sqrt( x_var_param(varidx,r-p1) ) ) -log(2*pi)/2 ...
                    - 0.5 * ( ( proposed - x_mean_param(varidx,r-p1) )./sqrt( x_var_param(varidx,r-p1) ) ).^2;
            %loglikeprop = loglikeprop + expit(sum(propTX_matrix .* beta_coeff(z(tmpidx,1),:),2));
            loglikeprop = loglikeprop + log(binopdf(Y(tmpidx,:),ones(sum(tmpidx),1), expit(sum(propTX_matrix .* beta_coeff(z(tmpidx,1),:),2)) ) + eps);
            loglikeprop = loglikeprop  -log( 0.1  ) -log(2*pi)/2 ...
                    - 0.5 * ( (X(tmpidx,r) -proposed )./ 0.1 ).^2;

            loglikeold = -log( sqrt( x_var_param(varidx,r-p1) ) ) -log(2*pi)/2 ...
                    - 0.5 * ( ( X(tmpidx,r) - x_mean_param(varidx,r-p1) )./sqrt( x_var_param(varidx,r-p1) ) ).^2;
%             loglikeold = loglikeold + expit(sum(TX_matrix(tmpidx,:) .* beta_coeff(z(tmpidx,1),:),2));
            loglikeold = loglikeold +  log(binopdf(Y(tmpidx,:),ones(sum(tmpidx),1), expit(sum(TX_matrix(tmpidx,:) .* beta_coeff(z(tmpidx,1),:),2)) ) + eps);
            loglikeold = loglikeold  -log( 0.1  ) -log(2*pi)/2 ...
                    - 0.5 * ( (proposed - X(tmpidx,r))./ 0.1 ).^2;
          
             indswitch = rand(size(loglikeold)) < exp(loglikeprop-loglikeold);
             X(tmpidx,r) = X(tmpidx,r) .* (1-indswitch) + indswitch .* proposed;
        end
    end
    
  % End update of missing covariates------------------------------
  
  % Now update variables that depend on updated covariates-----------------------  
  
  % Update design matrix
  % -------------------------------------------------------------------------
  X_matrix = [ones(size(T,1),1),X];
  TX = [T,X];
  TX_matrix = [ones(size(T,1),1),TX];
  end