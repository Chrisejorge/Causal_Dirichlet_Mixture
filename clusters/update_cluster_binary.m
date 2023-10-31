function [membership,betaY,betaT,x_pis,x_means,...
    x_vars] = update_cluster_binary(T,X,Y,p1,p2,z,beta_coeff,gamma_coeff,x_pi_param,x_mean_param,...
    x_var_param,alpha_theta,alpha_psi,alpha_omega,num_aux_param,beta_coeff0, diagbetacov0,...
    gamma_coeff0,diaggammacov0,c0, mu0, nu0, tau0, a0, b0)
% update clusters memberships

% Define design matrix
TX = [T,X];
TX_matrix = [ones(size(T,1),1),TX];
X_matrix = [ones(size(T,1),1),X];

% the indicators of clusters memberships for Y,T,X
zy = z(:,1);
zt = z(:,2);
zx = z(:,3);

% unique clusters memberships
uniqueZ = sortrows(unique(z,'rows'),1:size(z,2));
uniqueZyt= sortrows(unique(uniqueZ(:,1:2),'rows'),1:2);

% number of observations
n = size(X,1);

% For random generations of auxilary parameters
 m = num_aux_param;
 k = size(beta_coeff,1);
 beta_aux_table = mvnrnd(beta_coeff0,diagbetacov0, m * n );
%  beta_aux_table = mvnrnd(beta_coeff0,diagbetacov0, n );
 beta_table_idx = 1;
 gamma_aux_table = mvnrnd(gamma_coeff0,diaggammacov0, n * m * (1+ k) * 2);
%  gamma_aux_table = mvnrnd(gamma_coeff0,diaggammacov0,n );
 gamma_table_idx = 1;
 if p1 > 0
     xpipar_aux_table = betarnd(a0,b0,n * m * (1 + k + size(uniqueZyt,1)) * 2,p1);
%      xpipar_aux_table = betarnd(a0,b0,n,p1);
     xpipar_table_idx = 1;
 end
 if p2 > 0
     xvars_aux_table = sinvchi2rand(nu0,tau0*tau0,n * m * (1 + k + size(uniqueZyt,1)) * 2,p2);
%      xvars_aux_table = sinvchi2rand(nu0,tau0*tau0,n ,p2);
     xmupars_aux_table =  normrnd(mu0, sqrt(xvars_aux_table) / sqrt(c0));
     xvars_table_idx = 1;
 end

% to indicate whether the object is the only one
% for ii = 1:n
for ii = randperm(n)
    onlyT = 0;
    onlyY = 0;
    % check if the object is the only one in its cluster
    num_dummy = sum(zy == zy(ii) & zt == zt(ii) & zx == zx(ii));
    
    if num_dummy == 1
        % delete associated coefficients in the X cluster
        
        num_dummy2 = sum(zy == zy(ii) & zt == zt(ii));
        
        % check if the object is the only one in T cluster
        % delete T coef if only one in T cluster
        if num_dummy2 == 1
            onlyT = 1;
            
            num_dummy3 = sum(zy == zy(ii));
            
            if num_dummy3 == 1
                onlyY = 1;
                % check if the object is the only one in Y cluster
                % delete Y coef if only one in Y cluster
                beta_coeff(zy(ii),:) = [];
            end
            % delete T coef
            obj_t_idx = uniqueZyt(:,1)==zy(ii) & uniqueZyt(:,2)==zt(ii) ;
            gamma_coeff(obj_t_idx,:) = [];
        end % end num_dummy2
        % delete X coef
        obj_x_idx = uniqueZ(:,1)==zy(ii) & uniqueZ(:,2)==zt(ii) & uniqueZ(:,3)==zx(ii);
        if p1 > 0
        x_pi_param(obj_x_idx,:) = [];
        end
        if p2 > 0
        x_mean_param(obj_x_idx,:) = [];
        x_var_param(obj_x_idx,:) = [];
        end
        
        % relabel X cluster
        
        tmpidx = zy == zy(ii) & zt == zt(ii) & zx > zx(ii);
        zx(tmpidx) = zx(tmpidx)-1;
        
        tmpidx = uniqueZ(:,1) == zy(ii) & uniqueZ(:,2) == zt(ii) & uniqueZ(:,3) > zx(ii);
        uniqueZ(tmpidx,3) = uniqueZ(tmpidx,3)-1;
        
        % relabel T cluster
        if onlyT == 1
            
            tmpidx = zy == zy(ii) & zt > zt(ii);
            zt(tmpidx) = zt(tmpidx)-1;
            
            tmpidx = uniqueZ(:,1) == zy(ii) & uniqueZ(:,2) > zt(ii);
            uniqueZ(tmpidx,2) = uniqueZ(tmpidx,2)-1;
            
            tmpidx = uniqueZyt(:,1) == zy(ii) & uniqueZyt(:,2) > zt(ii);
            uniqueZyt(tmpidx,2) = uniqueZyt(tmpidx,2) - 1;
            
        end % end if onlyT
        
        % relabel Y cluster
        if onlyY == 1
            tmpidx = zy > zy(ii);
            zy(tmpidx) = zy(tmpidx)-1;
            
            tmpidx = uniqueZ(:,1) > zy(ii);
            uniqueZ(tmpidx,1) = uniqueZ(tmpidx,1)-1;
            
            tmpidx = uniqueZyt(:,1) > zy(ii);
            uniqueZyt(tmpidx,1) = uniqueZyt(tmpidx,1)-1;
        end % end if onlyY
        
        % remove object from the cluster uniquez
        uniqueZ(obj_x_idx,:) = [];
        
        % remove the object from the cluster uniqueZyt if needed
        if onlyT == 1
            uniqueZyt(obj_t_idx,:) = [];
        end
    end % end num_dummy
    
    % remove the object from zy,zt,zx
    zy(ii) = [];
    zt(ii) = [];
    zx(ii) = [];
    
    %% initilize parameters for calculation of probability
    % existing clusters plus auxiliary parameters;
    m = num_aux_param;
    k = size(beta_coeff,1);
%     assert(k == length(unique(zy)));
    if p2 > 0
        numTotalCluster = size(x_mean_param,1);
    else
        numTotalCluster = size(x_pi_param,1);
    end
    
    totalposs = numTotalCluster+m+k*m;
    kjs = zeros(k,1);
    for jj = 1:k
        kj = length( uniqueZyt(uniqueZyt(:,1) == jj,2));
        kjs(jj) = kj;
        totalposs = totalposs + kj*m;
    end
    
    probsWeights = zeros(totalposs,1);
    loglike = zeros(totalposs,1);
    count = 1;
    
    % init gamma_coeff_mat
%     maxNumTclusters = max(uniqueZyt(:,2));
%     gamma_coeff_mat = zeros(k+1,maxNumTclusters+1,size(X_matrix,2));
%     gcount = 1;
%     for jj = 1:k
%         tmpYidx = jj;
%         kj = kjs(jj);
%         for l = 1:kj
%             tmpXidx = l;
%             gamma_coeff_mat(tmpYidx,tmpXidx,:) = gamma_coeff(gcount,:);
%             gcount = gcount+1;
%         end
%     end % jj
    
    % to store Y,T,X indexs of existing clusters plus auxiliary parameters
    idxs_aux = ones(totalposs,3);
    
    %% calculate probs for clusters
            
    loglikeregy = 0;
    loglikeregt = 0;
    loglikeregx = 0;
    logprodx1 = 0;
    logprodx2 = 0;
    % probs, Y old T old X old
    loglikeregy_olds = logbernpdf(Y(ii),expit(sum(TX_matrix(ii,:) .* beta_coeff,2)));
    loglikeregt_olds = logbernpdf(T(ii),expit(sum(X_matrix(ii,:) .* gamma_coeff,2)));
    if p1 > 0
        logprodx1_olds = sum(logbernpdf(X(ii,1:p1),x_pi_param),2);
    end
    if p2 > 0
        logprodx2_olds = sum(lognormpdf(X(ii,p1+1:p1+p2), x_mean_param, sqrt( x_var_param ) ), 2 );
    end
    
    gcount = 1;
    x_param_count = 1;
    for jj = 1:k
        
        tmp_y_idx = jj;
        % get count of number of T clusters within jth Y cluster
        kj = kjs(jj);
        
        % get number of subjects within jth cluster
        njwoi = sum(zy == tmp_y_idx);
        
        % likelihood for each existing Y cluster
%         loglikeregy = log(binopdf(Y(ii),1,expit(TX_matrix(ii,:) * beta_coeff(jj,:)')) + eps);
%         loglikeregy = logbernpdf(Y(ii),expit(TX_matrix(ii,:) * beta_coeff(jj,:)'));
        loglikeregy = loglikeregy_olds(jj);
        for l = 1:kj
            tmp_t_idx = l;
            % get count of number of X clusters within jth Y cluster,lth T
            % cluster
            tmpidx = uniqueZ(:,1) == tmp_y_idx & uniqueZ(:,2) == tmp_t_idx;
            klj = length(uniqueZ(tmpidx,3));
            
%             tmpcoeff = reshape(gamma_coeff_mat(tmp_y_idx,tmp_t_idx,:),size(gamma_coeff0));
            
%             loglikeregt = log(binopdf(T(ii),1,expit(X_matrix(ii,:) * tmpcoeff')) + eps);
%             loglikeregt = logbernpdf(T(ii),expit(X_matrix(ii,:) * tmpcoeff'));
            loglikeregt = loglikeregt_olds(gcount);
            gcount = gcount+1;
            
            nljwoi = sum(zy == tmp_y_idx & zt == tmp_t_idx);
            
            for h = 1:klj
                tmp_x_idx = h;
                nhljwoi = sum(zy == tmp_y_idx & zt == tmp_t_idx  & zx == tmp_x_idx);
                                
                if p1 > 0
%                     logprodx1 = sum(log(binopdf(X(ii,1:p1),1,x_pi_param(x_param_count,:)) + eps));
%                     logprodx1 = sum(logbernpdf(X(ii,1:p1),x_pi_param(x_param_count,:)));
                    logprodx1 = logprodx1_olds(x_param_count);
                end
                if p2 > 0
%                     logprodx2 = sum(-log(sqrt(x_var_param(x_param_count,:))) -log(2*pi)/2 ...
%                         - 0.5 * ((X(ii,p1+1:p1+p2)-x_mean_param(x_param_count,:))./sqrt(x_var_param(x_param_count,:))).^2);
%                     logprodx2 = sum(lognormpdf(X(ii,p1+1:p1+p2), x_mean_param(x_param_count,:), sqrt( x_var_param( x_param_count,:) ) ) );
                    logprodx2 = logprodx2_olds(x_param_count);
                end
                
                loglikeregx = logprodx1 + logprodx2;
                
                x_param_count = x_param_count+1;
                
                idxs_aux(count,:) = [tmp_y_idx,tmp_t_idx,tmp_x_idx];
                
                probsWeights(count) = (njwoi*nljwoi*nhljwoi)/( ( n-1+alpha_theta) * (njwoi + alpha_psi) * (nljwoi + alpha_omega) );
                loglike(count) = loglikeregy + loglikeregt + loglikeregx;
                count = count + 1;
            end % for h
            
        end % for l
    end % for j
    
    % probs, Y old T old X new
    % set auxiliary parameters
    
    xcount = m*sum(kjs);
    
%     if p1 > 0
%         xpipar_aux = betarnd(a0,b0,xcount,p1);
%     end
%     if p2 > 0
%         xvars_aux = sinvchi2rand(nu0,tau0*tau0,xcount,p2);
%         xmupars_aux =  normrnd(mu0, sqrt(xvars_aux) / sqrt(c0));
%     end
    
    if p1 > 0
       if xpipar_table_idx + xcount-1 > size(xpipar_aux_table,1)
%            fprintf('xpipar fantasy \n');
           xpipar_aux_table = betarnd(a0,b0,n * m * k * k * 3,p1);
           xpipar_table_idx = 1;
       end
       xpipar_aux = xpipar_aux_table(xpipar_table_idx:xpipar_table_idx+xcount-1,:);
       xpipar_table_idx = xpipar_table_idx + xcount;
    end
    if p2 > 0
        if xvars_table_idx + xcount-1 > size(xvars_aux_table,1)
%            fprintf('xvars fantasy \n');
            xvars_aux_table = sinvchi2rand(nu0,tau0*tau0,n * m * (1 + k + size(uniqueZyt,1)) * 2,p2);
            xmupars_aux_table =  normrnd(mu0, sqrt(xvars_aux_table) / sqrt(c0));
            xvars_table_idx = 1;
       end
       xvars_aux = xvars_aux_table(xvars_table_idx:xvars_table_idx+xcount-1,:);
       xmupars_aux = xmupars_aux_table(xvars_table_idx:xvars_table_idx+xcount-1,:);
       xvars_table_idx = xvars_table_idx + xcount;
    end
    
    if p1 > 0
        logprodx1_new= sum(logbernpdf(X(ii,1:p1),xpipar_aux),2);
    end
    
    if p2 > 0
        logprodx2_new = sum(lognormpdf(X(ii,p1+1:p1+p2), xmupars_aux, sqrt( xvars_aux ) ), 2 );
    end
    
    gcount = 1;
    xpar_aux_count = 1;
    newXidx = count;
    for jj = 1:k
        
        tmp_y_idx = jj;
        % get count of number of T clusters within jth Y cluster
        kj = kjs(jj);
        
        % get number of subjects within jth cluster
        njwoi = sum(zy == tmp_y_idx);
        
        % likelihood for each existing Y cluster
%         loglikeregy = log(binopdf(Y(ii),1,expit(TX_matrix(ii,:) * beta_coeff(jj,:)')) + eps);
%         loglikeregy = logbernpdf(Y(ii),expit(TX_matrix(ii,:) * beta_coeff(jj,:)'));
        loglikeregy = loglikeregy_olds(jj);
        for l = 1:kj
            tmp_t_idx = l;
            
%              tmpcoeff = reshape(gamma_coeff_mat(tmp_y_idx,tmp_t_idx,:),size(gamma_coeff0));
% % %             loglikeregt = log(binopdf(T(ii),1,expit(X_matrix(ii,:) * tmpcoeff')) + eps);
%              loglikeregt = logbernpdf(T(ii),expit(X_matrix(ii,:) * tmpcoeff'));
            
            loglikeregt = loglikeregt_olds(gcount);
            gcount = gcount+1;
            
            nljwoi = sum(zy == tmp_y_idx & zt == tmp_t_idx);
            
            for h = 1:m
                             
                if p1 > 0
%                     logprodx1 = sum(log(binopdf(X(ii,1:p1),1,xpipar_aux(xpar_aux_count,:)) + eps));
%                     logprodx1 = sum(logbernpdf(X(ii,1:p1),xpipar_aux(xpar_aux_count,:)));
                    logprodx1 = logprodx1_new(xpar_aux_count);
                end
                
                if p2 > 0
%                     logprodx2 = sum(-log(sqrt(xvars_aux(xpar_aux_count,:))) -log(2*pi)/2 ...
%                         - 0.5 * ((X(ii,p1+1:p1+p2)-xmupars_aux(xpar_aux_count,:))./sqrt(xvars_aux(xpar_aux_count,:))).^2);
%                     logprodx2 = sum(lognormpdf(X(ii,p1+1:p1+p2), xmupars_aux(xpar_aux_count,:), sqrt( xvars_aux(xpar_aux_count,:)) ) );
                    logprodx2 = logprodx2_new(xpar_aux_count);
                end
                
                loglikeregx = logprodx1 + logprodx2;
                
                xpar_aux_count = xpar_aux_count + 1;
                
                x_idx_aux = length(uniqueZ(uniqueZ(:,1) == tmp_y_idx & uniqueZ(:,2) == tmp_t_idx,3) ) + 1;
                idxs_aux(count,:) = [tmp_y_idx,tmp_t_idx,x_idx_aux];
                
                probsWeights(count) = (njwoi*nljwoi*alpha_omega/m)/( ( n-1+alpha_theta) * (njwoi + alpha_psi) * (nljwoi + alpha_omega) );
                loglike(count) = loglikeregy + loglikeregt + loglikeregx;
                count = count + 1;
            end % for h
            
        end % for l
    end %  for j
    
    % probs, Y old T new X new
    % set auxiliary parameters
    
%     tcount = k*m;
%     gamma_coeff_aux = mvnrnd(gamma_coeff0,diaggammacov0,tcount);
    
     tcount = k*m;
     if gamma_table_idx + tcount-1 > size(gamma_aux_table,1)
%         fprintf('gamma fantasy ');
        gamma_aux_table = mvnrnd(gamma_coeff0,diaggammacov0, n * m * k * 3);
        gamma_table_idx = 1;
    end
    gamma_coeff_aux = gamma_aux_table(gamma_table_idx:gamma_table_idx+tcount-1,:);
    gamma_table_idx = gamma_table_idx + tcount;
    
    loglikeregt_new = logbernpdf(T(ii),expit(sum(X_matrix(ii,:) .* gamma_coeff_aux,2)));
    
%     if p1 > 0
%         t_xpipar_aux = betarnd(a0,b0,tcount,p1);
%     end
%     if p2 > 0
%         t_xvars_aux = sinvchi2rand(nu0,tau0*tau0,tcount,p2);
%         t_xmupars_aux = normrnd(mu0, sqrt(t_xvars_aux)/sqrt(c0));
%     end
    
     if p1 > 0
       if xpipar_table_idx + tcount-1 > size(xpipar_aux_table,1)
%            fprintf('xpipar fantasy \n');
           xpipar_aux_table = betarnd(a0,b0,n * m * k * k * 3,p1);
           xpipar_table_idx = 1;
       end
       t_xpipar_aux = xpipar_aux_table(xpipar_table_idx:xpipar_table_idx+tcount-1,:);
       xpipar_table_idx = xpipar_table_idx + tcount;
     end
    
     if p2 > 0
        if xvars_table_idx + tcount-1 > size(xvars_aux_table,1)
%            fprintf('xvars fantasy \n');
            xvars_aux_table = sinvchi2rand(nu0,tau0*tau0,n * m * (1 + k + size(uniqueZyt,1)) * 2,p2);
            xmupars_aux_table =  normrnd(mu0, sqrt(xvars_aux_table) / sqrt(c0));
            xvars_table_idx = 1;
       end
       t_xvars_aux = xvars_aux_table(xvars_table_idx:xvars_table_idx+tcount-1,:);
       t_xmupars_aux = xmupars_aux_table(xvars_table_idx:xvars_table_idx+tcount-1,:);
       xvars_table_idx = xvars_table_idx + tcount;
    end
    
    if p1 > 0
        logprodx1_new= sum(logbernpdf(X(ii,1:p1),t_xpipar_aux),2);
    end
    
    if p2 > 0
        logprodx2_new = sum(lognormpdf(X(ii,p1+1:p1+p2), t_xmupars_aux, sqrt( t_xvars_aux ) ), 2 );
    end
    
    newTidx = count;
    t_aux_count = 1;
    for jj = 1:k
        
        tmp_y_idx = jj;
        
        % get number of subjects within jth cluster
        njwoi = sum(zy == tmp_y_idx);
        
        % likelihood for each existing Y cluster
%         loglikeregy = log(binopdf(Y(ii),1,expit(TX_matrix(ii,:) * beta_coeff(jj,:)')) + eps);
%         loglikeregy = logbernpdf(Y(ii),expit(TX_matrix(ii,:) * beta_coeff(jj,:)'));
        loglikeregy = loglikeregy_olds(jj);
        for l = 1:m
            
%             loglikeregt = log(binopdf(T(ii),1,expit(X_matrix(ii,:) * gamma_coeff_aux(t_aux_count,:)')) + eps);
%             loglikeregt = logbernpdf(T(ii),expit(X_matrix(ii,:) * gamma_coeff_aux(t_aux_count,:)'));
            loglikeregt = loglikeregt_new(t_aux_count);
                     
            if p1 > 0
%                 logprodx1 = sum(log(binopdf(X(ii,1:p1),1,t_xpipar_aux(t_aux_count,:)) + eps));
%                 logprodx1 = sum(logbernpdf(X(ii,1:p1),t_xpipar_aux(t_aux_count,:)));
                logprodx1 = logprodx1_new(t_aux_count);
            end
            
            if p2 > 0
%                 logprodx2 = sum(-log(sqrt(t_xvars_aux(t_aux_count,:))) -log(2*pi)/2 ...
%                     - 0.5 * ((X(ii,p1+1:p1+p2)-t_xmupars_aux(t_aux_count,:))./sqrt(t_xvars_aux(t_aux_count,:))).^2);
%                logprodx2 = sum(lognormpdf(X(ii,p1+1:p1+p2), t_xmupars_aux(t_aux_count,:), sqrt( t_xvars_aux(t_aux_count,:)) ) );
               logprodx2 = logprodx2_new(t_aux_count);
            end
            
            loglikeregx = logprodx1 + logprodx2;
            
            t_aux_count = t_aux_count + 1;
            
            idxs_aux(count,:) = [tmp_y_idx,kjs(jj)+1,1];
            
            probsWeights(count) = (njwoi*alpha_psi/m)/( ( n-1+alpha_theta) * (njwoi + alpha_psi) );
            loglike(count) = loglikeregy + loglikeregt + loglikeregx;
            count = count + 1;
        end % for l
        
    end % for j
    
    % probs, Y new (T,X new)
    % set auxiliary parameters
    
%     beta_coeff_aux = mvnrnd(beta_coeff0,diagbetacov0,m);
      if beta_table_idx + m-1 > size(beta_aux_table,1)
%            fprintf('beta fantasy');
           beta_aux_table = mvnrnd(beta_coeff0,diagbetacov0, m * n );
           beta_table_idx = 1;
      end
      beta_coeff_aux = beta_aux_table(beta_table_idx:beta_table_idx+m-1,:);
      beta_table_idx = beta_table_idx + m;
      
%      y_gamma_coeff_aux = mvnrnd(gamma_coeff0,diaggammacov0,m);

    if gamma_table_idx + m-1 > size(gamma_aux_table,1)
%         fprintf('gamma fantasy');
        gamma_aux_table = mvnrnd(gamma_coeff0,diaggammacov0, n * m * k * 3);
        gamma_table_idx = 1;
    end
    y_gamma_coeff_aux = gamma_aux_table(gamma_table_idx:gamma_table_idx+m-1,:);
    gamma_table_idx = gamma_table_idx + m;
    
%     if p1 > 0
%         y_xpipar_aux = betarnd(a0,b0,m,p1);
%     end
%     if p2 > 0
%         y_xvars_aux = sinvchi2rand(nu0,tau0*tau0,m,p2);
%         y_xmupars_aux =  normrnd(mu0, sqrt(y_xvars_aux)/sqrt(c0));
%     end
    
     if p1 > 0
       if xpipar_table_idx + m-1 > size(xpipar_aux_table,1)
%            fprintf('xpipar fantasy \n');
           xpipar_aux_table = betarnd(a0,b0,n * m * k * k * 3,p1);
           xpipar_table_idx = 1;
       end
       y_xpipar_aux = xpipar_aux_table(xpipar_table_idx:xpipar_table_idx+m-1,:);
       xpipar_table_idx = xpipar_table_idx + m;
     end
    
     if p2 > 0
        if xvars_table_idx + m-1 > size(xvars_aux_table,1)
%            fprintf('xvars fantasy \n');
            xvars_aux_table = sinvchi2rand(nu0,tau0*tau0,n * m * (1 + k + size(uniqueZyt,1)) * 2,p2);
            xmupars_aux_table =  normrnd(mu0, sqrt(xvars_aux_table) / sqrt(c0));
            xvars_table_idx = 1;
       end
       y_xvars_aux = xvars_aux_table(xvars_table_idx:xvars_table_idx+m-1,:);
       y_xmupars_aux = xmupars_aux_table(xvars_table_idx:xvars_table_idx+m-1,:);
       xvars_table_idx = xvars_table_idx + m;
    end
    
    loglikeregy_new = logbernpdf(Y(ii),expit(sum(TX_matrix(ii,:) .* beta_coeff_aux,2)));
    loglikeregt_new = logbernpdf(T(ii),expit(sum(X_matrix(ii,:) .* y_gamma_coeff_aux,2)));
    
     if p1 > 0
        logprodx1_new= sum(logbernpdf(X(ii,1:p1),y_xpipar_aux),2);
    end
    
    if p2 > 0
        logprodx2_new = sum(lognormpdf(X(ii,p1+1:p1+p2), y_xmupars_aux, sqrt( y_xvars_aux ) ), 2 );
    end
    
    newYidx = count;
    for jj = 1:m
        
        % likelihood for each new Y cluster
%         loglikeregy = log(binopdf(Y(ii),1,expit(TX_matrix(ii,:) * beta_coeff_aux(jj,:)')) + eps);
%         loglikeregy = logbernpdf(Y(ii),expit(TX_matrix(ii,:) * beta_coeff_aux(jj,:)'));
        loglikeregy = loglikeregy_new(jj);
        
%         loglikeregt = log( binopdf(T(ii),1,expit(X_matrix(ii,:) * y_gamma_coeff_aux(jj,:)')) + eps);
%         loglikeregt = logbernpdf(T(ii),expit(X_matrix(ii,:) * y_gamma_coeff_aux(jj,:)'));
        loglikeregt = loglikeregt_new(jj);
        
        if p1 > 0
%             logprodx1 = sum(log(binopdf(X(ii,1:p1),1,y_xpipar_aux(jj,:)) + eps));
%             logprodx1 = sum(logbernpdf(X(ii,1:p1),y_xpipar_aux(jj,:)));
            logprodx1 = logprodx1_new(jj);
        end
        
        if p2 > 0
%             logprodx2 = sum(-log(sqrt(y_xvars_aux(jj,:))) -log(2*pi)/2 ...
%                 - 0.5 * ((X(ii,p1+1:p1+p2)-y_xmupars_aux(jj,:))./sqrt(y_xvars_aux(jj,:))).^2);
%             logprodx2 = sum(lognormpdf(X(ii,p1+1:p1+p2), y_xmupars_aux(jj,:), sqrt( y_xvars_aux(jj,:) ) ) );
            logprodx2 = logprodx2_new(jj);
        end
        
        loglikeregx = logprodx1 + logprodx2;
        
        idxs_aux(count,:) = [k+1,1,1];
        
        probsWeights(count) = (alpha_theta/m)/ ( n-1+alpha_theta);
        loglike(count) = loglikeregy + loglikeregt + loglikeregx;
        count = count + 1;
    end % for j
    
    % Calculation of probs
    probs = probsWeights.* exp(loglike-max(loglike));
%     probs = probsWeights;
    % select cluster
    newCluster = sum (rand > cumsum(probs/sum(probs)))  + 1;
    %clear probs;
    % update parameters of clusters
    if newCluster < newXidx
        % old old old
        % fprintf( [ 'old old old\n' ] );
        tmpYidx = idxs_aux(newCluster,1);
        tmpTidx = idxs_aux(newCluster,2);
        tmpXidx = idxs_aux(newCluster,3);
        zy = [zy(1:ii-1);tmpYidx;zy(ii:end)];
        zt = [zt(1:ii-1);tmpTidx;zt(ii:end)];
        zx = [zx(1:ii-1);tmpXidx;zx(ii:end)];
        %         assert(idxs_aux(newCluster,1) == uniqueZ(newCluster,1));
        %         assert(idxs_aux(newCluster,2) == uniqueZ(newCluster,2));
        %         assert(idxs_aux(newCluster,3) == uniqueZ(newCluster,3));
    else
        if  newCluster >=  newYidx
            % new new new
            
            tmpYidx = idxs_aux(newCluster,1);
            tmpTidx = idxs_aux(newCluster,2);
            tmpXidx = idxs_aux(newCluster,3);
            
            zy = [zy(1:ii-1);tmpYidx;zy(ii:end)];
            zt = [zt(1:ii-1);tmpTidx;zt(ii:end)];
            zx = [zx(1:ii-1);tmpXidx;zx(ii:end)];
            
            %             assert(idxs_aux(newCluster,1) == k+1);
            %             assert(idxs_aux(newCluster,2) == 1);
            %             assert(idxs_aux(newCluster,3) == 1);
            
            uniqueZ = [uniqueZ;[zy(ii),zt(ii),zx(ii)]];
            uniqueZyt = [uniqueZyt;[zy(ii),zt(ii)]];
            
            tmpidx_aux = newCluster-newYidx+1;
            beta_coeff = [beta_coeff;beta_coeff_aux(tmpidx_aux,:)];
            gamma_coeff = [gamma_coeff; y_gamma_coeff_aux(tmpidx_aux,:)];
            
            if p1 > 0
                x_pi_param = [x_pi_param; y_xpipar_aux(tmpidx_aux,:)];
            end
            
            if p2 > 0
                x_mean_param = [x_mean_param; y_xmupars_aux(tmpidx_aux,:)];
                x_var_param = [x_var_param; y_xvars_aux(tmpidx_aux,:)];
            end
            
        else
            if newCluster >= newTidx
                %old new new
                %                  fprintf( [ 'old new new\n' ] );
                
                %                 tmpYidx = ceil((newCluster-newTidx+1)/m);
                % tmpTidx = length(uniqueZyt(uniqueZyt(:,1) == tmpYidx,2))+1;
                %                 tmpXidx = 1;
                
                tmpYidx = idxs_aux(newCluster,1);
                zy = [zy(1:ii-1);tmpYidx;zy(ii:end)];
                %                 assert(tmpYidx == ceil((newCluster-newTidx+1)/m));
                
                
                tmpTidx = idxs_aux(newCluster,2);
                zt = [zt(1:ii-1);tmpTidx;zt(ii:end)];
                %                 assert(tmpTidx == length(uniqueZyt(uniqueZyt(:,1) == tmpYidx,2))+1);
                
                tmpXidx = idxs_aux(newCluster,3);
                zx = [zx(1:ii-1);tmpXidx;zx(ii:end)];
                %                 assert(tmpXidx == 1);
                
                tmpidx_aux = newCluster-newTidx+1;
                tmpidx = sum(uniqueZ(:,1) <= tmpYidx)+1;
                uniqueZ = [uniqueZ(1:tmpidx-1,:); [tmpYidx,tmpTidx,tmpXidx];uniqueZ(tmpidx:end,:)];
                if p1 > 0
                    x_pi_param = [ x_pi_param(1:tmpidx-1,:);t_xpipar_aux(tmpidx_aux,:);x_pi_param(tmpidx:end,:) ];
                end
                
                if p2 > 0
                    x_var_param = [ x_var_param(1:tmpidx-1,:);t_xvars_aux(tmpidx_aux,:);x_var_param(tmpidx:end,:) ];
                    x_mean_param = [ x_mean_param(1:tmpidx-1,:);t_xmupars_aux(tmpidx_aux,:);x_mean_param(tmpidx:end,:) ];
                end
                
                tmpidx = length(find(uniqueZyt(:,1) <= tmpYidx))+1;
                uniqueZyt = [uniqueZyt(1:tmpidx-1,:); [tmpYidx,tmpTidx];uniqueZyt(tmpidx:end,:) ];
                gamma_coeff = [gamma_coeff(1:tmpidx-1,:);gamma_coeff_aux(tmpidx_aux,:); gamma_coeff(tmpidx:end,:)];
                
                
                
            elseif newCluster >= newXidx
                % old old new
                %                  fprintf( [ 'old old new\n' ] );
                
                % get count of number of T clusters within jth Y cluster
                %                  yidxs = crosstab(uniqueZyt(:,1)) * m;
                %                  yidxs = cumsum(yidxs);
                
                %tmpYidx = sum (newCluster-newXidx+1 > yidxs)  + 1;
                %tmpTidx = ceil((newCluster-newXidx+1-tmp)/m);
                % tmpXidx = size(uniqueZ(uniqueZ(:,1) == tmpYidx & uniqueZ(:,2) == tmpTidx),1) + 1;
                
                
                tmpYidx = idxs_aux(newCluster,1);
                zy = [zy(1:ii-1);tmpYidx;zy(ii:end)];
                %                 assert(tmpYidx == sum (newCluster-newXidx+1 > yidxs)  + 1);
                
                %                 if tmpYidx < 2
                %                     tmp = 0;
                %                 else
                %                     tmp = yidxs(tmpYidx-1);
                %                 end
                
                tmpTidx = idxs_aux(newCluster,2);
                zt = [zt(1:ii-1);tmpTidx;zt(ii:end)];
                %                 assert(tmpTidx == ceil((newCluster-newXidx+1-tmp)/m));
                
                tmpXidx = idxs_aux(newCluster,3);
                zx = [zx(1:ii-1);tmpXidx;zx(ii:end)];
                %                 assert(tmpXidx == size(uniqueZ(uniqueZ(:,1) == tmpYidx & uniqueZ(:,2) == tmpTidx),1) + 1);
                
                
                tmpidx_aux = newCluster-newXidx+1;
                tmpidx = sum(uniqueZ(:,1) <= tmpYidx) - sum(uniqueZ(uniqueZ(:,1) == tmpYidx,2) > tmpTidx)+1;
                %                 assert(sum(uniqueZ(:,1) <= tmpYidx) - sum(uniqueZ(uniqueZ(:,1) == tmpYidx,2) > tmpTidx) == sum(uniqueZ(:,1) < tmpYidx) + sum(uniqueZ(uniqueZ(:,1) == tmpYidx,2) <= tmpTidx) );
                uniqueZ = [uniqueZ(1:tmpidx-1,:); [tmpYidx,tmpTidx,tmpXidx]; uniqueZ(tmpidx:end,:)];
                
                if p1 > 0
                    x_pi_param = [x_pi_param(1:tmpidx-1,:); xpipar_aux(tmpidx_aux,:);x_pi_param(tmpidx:end,:)];
                end
                
                if p2 > 0
                    x_mean_param = [x_mean_param(1:tmpidx-1,:); xmupars_aux(tmpidx_aux,:);x_mean_param(tmpidx:end,:)];
                    x_var_param = [x_var_param(1:tmpidx-1,:); xvars_aux(tmpidx_aux,:);x_var_param(tmpidx:end,:)];
                end
                
            end % old new new
        end % new new new
        
    end %  old old old
    
end % end cluster loop

membership = [zy,zt,zx];
betaY = beta_coeff;
betaT = gamma_coeff;
x_pis = x_pi_param;
x_means = x_mean_param;
x_vars = x_var_param;

% end
