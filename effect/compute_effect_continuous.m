function new_psi_rd = compute_effect_continuous(X,n,p1,p2,z,k,kjs,kljs,njs,nljs,nhljs,beta_coeff,gamma_coeff,x_pi_param,x_mean_param,x_var_param,...
    alpha_theta,alpha_psi,alpha_omega,beta_coeff0,diagbetacov0,gamma_coeff0,diaggammacov0,c0, mu0, nu0, tau0, a0, b0,M,beta_coeff_prior_mc,gamma_coeff_prior_mc)
% Calculation of causal effects

gamma_coeff_mat = get_gamma_coeff_mat(z,k,kjs,gamma_coeff);

if p1 > 0
    x_pi_param_mat = get_pi_param_mat(z,k,kjs,kljs,x_pi_param);
end

if p2 > 0
    x_mean_param_mat = get_mean_param_mat(z,k,kjs,kljs,x_mean_param);
    x_var_param_mat = get_var_param_mat(z,k,kjs,kljs,x_var_param);
end

% Initialize matrices and vectors
        beta_coeff_mc = zeros(k+1,size(beta_coeff0,2));
        beta_coeff_mc(1:k,:) = beta_coeff;
              
        z_mc = ones(M,3);
        X_mc = NaN(M,size(X,2));
        
        % Assign M observations one of the current y-clusters or to a new
        % y-cluster using draws from a multinomial distribution with probabilities 'prob_ycluster'
        
        pdenom = n + alpha_theta;
        prob_ycluster = [njs;alpha_theta]./pdenom;
        
        z_mc(:,1) = randsample(1:length(prob_ycluster),M,true,prob_ycluster);
        
        % Store indices of non-empty y-clusters (some y-clusters could be empty and there
        % could be a new cluster)
        yclusters_mc = sort(unique(z_mc(:, 1)));
        
        % Store number of y clusters including empty ones and possibly a new one
        k_mc = max(yclusters_mc);
        
        % Store number of observations in each y-cluster (should add to a total of M)
        Njs_mc = zeros(k_mc,0);
        
        % If a new y-cluster was opened, draw parameters for this cluster from priors.
        % new new new
        if k_mc > k
            gamma_coeff_mat(k_mc,1,:) = mvnrnd(gamma_coeff0,diaggammacov0,1);
                      
            if p1 > 0
                x_pi_param_mat(k_mc,1,1,:) = betarnd(a0,b0,1,p1);
            end
            
            if p2 > 0
                x_var_param_mat(k_mc,1,1,:)  = sinvchi2rand(nu0,tau0*tau0,1,p2);
                tmp_x_var_param = reshape(x_var_param_mat(k_mc,1,1,:),1,p2);
                x_mean_param_mat(k_mc,1,1,:) = normrnd(mu0,sqrt(tmp_x_var_param/c0));
            end
            
            beta_coeff_mc(k_mc,:) = mvnrnd(beta_coeff0,diagbetacov0,1);
        end
        
         % Draw t-clusters within each y-cluster
        % old new new
        for jj = yclusters_mc'
            % Number of observations (out of M total) in each y-cluster
            Njs_mc(jj) = sum(z_mc(:,1) == jj);
            if jj < k+1
                kj = kjs(jj);
                
                ptdenom = njs(jj) + alpha_psi;
                ptcluster = [nljs(jj,1:kj)';alpha_psi]./ptdenom;
                
                z_mc(z_mc(:,1) == jj,2) = randsample(1:length(ptcluster),Njs_mc(jj),true,ptcluster);
                % If a new t-cluster was opened, draw parameter from prior
                % old new new
                if max(z_mc(z_mc(:,1) == jj,2)) == kj+1
                    tmpidx = kj + 1;
                    gamma_coeff_mat(jj,tmpidx,:) = mvnrnd(gamma_coeff0,diaggammacov0,1);
                    
                    % parameter for x-cluster
                    if p1 > 0                      
                        x_pi_param_mat(jj,tmpidx,1,:) = betarnd(a0,b0,1,p1);
                    end
                    
                    if p2 > 0                     
                        x_var_param_mat(jj,tmpidx,1,:)  = sinvchi2rand(nu0,tau0*tau0,1,p2);
                        tmp_x_var_param = reshape(x_var_param_mat(jj,tmpidx,1,:),1,p2);
                        x_mean_param_mat(jj,tmpidx,1,:) = normrnd(mu0,sqrt(tmp_x_var_param/c0));
                    end
                end
            end
            % z_mc of t set to one when jj == k+1, we have done that
        end
        
         % Store number of observations in each t-cluster (should add to a total of M)
        Nljs_mc = zeros(k_mc,max(z_mc(:,2)));
        % Draw x-clusters within each t-cluster
        % old old new
        for jj = yclusters_mc'
            if jj < k+1
                kj = kjs(jj);
                % Store indices of non-empty t-clusters (some t-clusters could be empty and there
                % could be a new cluster)
                tcluster_mc = sort(unique(z_mc(z_mc(:,1) == jj,2)));
                for l = tcluster_mc'
                    tmpidx = z_mc(:,1) == jj & z_mc(:,2) == l;
                    Nljs_mc(jj,l) = sum(tmpidx);
                    
                    if l < kj+1
                        klj = kljs(jj,l);
                        
                        pxdenom = nljs(jj,l) + alpha_omega;
                        pxcluster = [reshape(nhljs(jj,l,1:klj),klj,1);alpha_omega]./pxdenom;
                        
                        z_mc(tmpidx,3) = randsample(1:length(pxcluster),Nljs_mc(jj,l),true,pxcluster);
                        
                        % if a new x-cluster was opened, draw parameter from prior
                        if max(z_mc(tmpidx,3)) == klj+1
                            tmp = klj+1;
                            if p1 > 0                               
                                x_pi_param_mat(jj,l,tmp,:) = betarnd(a0,b0,1,p1);
                            end
                            
                            if p2 > 0                     
                                x_var_param_mat(jj,l,tmp,:)   = sinvchi2rand(nu0,tau0*tau0,1,p2);
                                tmp_x_var_param = reshape(x_var_param_mat(jj,l,tmp,:),1,p2);
                                x_mean_param_mat(jj,l,tmp,:) = normrnd(mu0,sqrt(tmp_x_var_param/c0));
                            end
                            
                        end
                        
                    end
                end % l
            else
                tcluster_mc = sort(unique(z_mc(z_mc(:,1) == jj,2)));
                for l = tcluster_mc'
                    tmpidx = z_mc(:,1) == jj & z_mc(:,2) == l;
                    Nljs_mc(jj,l) = sum(tmpidx);
                end
            end %  if jj
            
        end % for jj
        
        assert(M == sum(sum(Nljs_mc)));
        
        Nhljs_mc = zeros(k_mc,max(z_mc(:,2)),max(z_mc(:,3)));
        % Draw covariates for each of the M observations-------------------
        for jj = yclusters_mc'
            tcluster_mc = sort(unique(z_mc(z_mc(:,1) == jj,2)));
            for l = tcluster_mc'
                tmpidx1 = z_mc(:,1) == jj & z_mc(:,2) == l;
                xcluster_mc = sort(unique(z_mc( tmpidx1,3) ) );
                for h = xcluster_mc'
                    tmpidx = tmpidx1 & z_mc(:,3) == h;
                    Nhljs_mc(jj,l,h) = sum(tmpidx);
                    
                    if p1 > 0                     
                        tmp_pi_param = reshape(x_pi_param_mat(jj,l,h,:),1,p1);
                        X_mc(tmpidx,1:p1) =  binornd(1,repmat(tmp_pi_param,Nhljs_mc(jj,l,h),1));
                    end
                    
                    if p2 > 0                      
                        tmp_x_var_param = reshape(x_var_param_mat(jj,l,h,:),1,p2);
                        tmp_x_mean_param = reshape(x_mean_param_mat(jj,l,h,:),1,p2);
                        X_mc(tmpidx,p1+1:p1+p2) =  normrnd(repmat(tmp_x_mean_param,Nhljs_mc(jj,l,h),1),repmat(sqrt(tmp_x_var_param),Nhljs_mc(jj,l,h),1));
                    end
                end
            end % for l
        end % for jj       
        
         % Compute E(Y|T=0, X=x,theta,psi,omega,z) and E(Y|T=1, X=x,theta,psi,omega,z) using M new
        % draws for covariates.-------------------------------------------------------
        
        f0Xt1_mat_mc = zeros(size(X_mc));
        f0Xt0_mat_mc = zeros(size(X_mc));
        % F_{0}(x) for T = 1 and T = 0
        f0Xt1_mc = zeros(size(X_mc,1),1);
        f0Xt0_mc = zeros(size(X_mc,1),1);
        
        % compute F_{0}(x)
        
        beta_a0b0  = beta(a0,b0);
        beta_draw = beta(a0 + X_mc(:,1:p1), b0 - X_mc(:,1:p1)+1)./beta_a0b0;
        f0Xt1_mat_mc(:,1:p1) = beta_draw;
        f0Xt0_mat_mc(:,1:p1) = beta_draw;
                
        cn = c0 + 1;
        nun = nu0 + 1;
        num_mc =  gamma( nun / 2 )*( nu0 * tau0 * tau0 ) ^ ( nu0 / 2 );
        tmpx = X_mc(:,p1+1:p1+p2);
        taun2_mc = (1/nun) * (nu0 * tau0 * tau0 + (c0/cn) .* (mu0 - tmpx).^2);
        denom = gamma(nu0/2) .* (nun .* taun2_mc).^(nun/2);
        frac_mc = num_mc * sqrt(c0/(cn * pi)) ./ denom;
        f0Xt1_mat_mc(:,p1+1:p1+p2) = frac_mc;
        f0Xt0_mat_mc(:,p1+1:p1+p2) = frac_mc;
        
        % Take product (covariates are assumed to be locally independent).
        % Result is vector of size M
        f0Xt1_mc(:,1) = prod(f0Xt1_mat_mc,2);
        f0Xt0_mc(:,1) = prod(f0Xt0_mat_mc,2);
        
        % compute F_{0}(t|x)
        X_matrix_mc = [ones(size(X_mc,1),1),X_mc];
        
%         f0t1_mc = mean(binopdf(ones(size(X_matrix_mc,1),size(gamma_coeff_prior_mc,1)),1,expit(X_matrix_mc * gamma_coeff_prior_mc')),2);
%         f0t0_mc = mean(binopdf(zeros(size(X_matrix_mc,1),size(gamma_coeff_prior_mc,1)),1,expit(X_matrix_mc * gamma_coeff_prior_mc')),2);
        
        tmp = expit(X_matrix_mc * gamma_coeff_prior_mc');
        f0t1_mc = mean(tmp,2);
        f0t0_mc = mean(1-tmp,2);
        
        % E_{0}[y|t,x]
        % Set treatment to 0 and 1, and add column of 1's.
        TXt1_matrix_mc = [ones(size(X_mc,1),1),ones(size(X_mc,1),1),X_mc];
        TXt0_matrix_mc = [ones(size(X_mc,1),1),zeros(size(X_mc,1),1),X_mc];
        
        e0yt1_mc = mean(TXt1_matrix_mc * beta_coeff_prior_mc',2);
        e0yt0_mc = mean(TXt0_matrix_mc * beta_coeff_prior_mc',2);
        % Calculate causal effects --------------------------------------
        % Calculate E[Y^t] using MC integration to integrate E(Y|T=t,X=x,theta^{*},psi^{*},omega^{*},z)
        % over covariate distribution
        partt1_mc = zeros(size(X_mc,1),1);
        partt0_mc = zeros(size(X_mc,1),1);
        wt1_mc = zeros(size(X_mc,1),k_mc);
        wt0_mc = zeros(size(X_mc,1),k_mc);
        for jj = yclusters_mc'
            % (Use yclusters.mc rather than 1:k.mc in case some y-clusters are empty)
            tcluster_mc = sort(unique(z_mc(z_mc(:,1) == jj,2)));
            nj_mc = Njs_mc(jj);
            tmpmax = max(tcluster_mc);
            wwt1_mc = zeros(size(X_mc,1),tmpmax);
            wwt0_mc = zeros(size(X_mc,1),tmpmax);
            sumt1_l_mc = zeros(size(X_mc,1),1);
            sumt0_l_mc = zeros(size(X_mc,1),1);
            for l = tcluster_mc'
                tmpidx1 = z_mc(:,1) == jj & z_mc(:,2) == l;
                xcluster_mc = sort(unique(z_mc( tmpidx1,3) ) );
                nlj_mc = Nljs_mc(jj,l);
                
                sum_h_mc = zeros(size(X_mc,1),1);
                for h   = xcluster_mc'
                    % F(x|w^{*})%
                    logprodx1_mc = 0;
                    if p1 > 0
                        %
                        tmp_pi_param = reshape(x_pi_param_mat(jj,l,h,:),1,p1);
%                         logprodx1_mc =  sum( log(binopdf(X_mc(:,1:p1),ones(size(X_mc,1),p1),repmat(tmp_pi_param,size(X_mc,1),1) ) + eps),2 ) ;
                        logprodx1_mc = sum(logbernpdf(X_mc(:,1:p1),tmp_pi_param),2);
                    end
                    
                    logprodx2_mc = 0;
                    if p2 > 0
                        tmp_x_var_param = reshape(x_var_param_mat(jj,l,h,:),1,p2);
                        tmp_x_mean_param = reshape(x_mean_param_mat(jj,l,h,:),1,p2);
                        logprodx2_mc = sum(-log(sqrt(tmp_x_var_param)) -log(2*pi)/2 ...
                            - 0.5 * ((X_mc(:,p1+1:p1+p2)-tmp_x_mean_param)./sqrt(tmp_x_var_param)).^2,2);
                        
                    end
                    nhjl_mc = Nhljs_mc(jj,l,h);
                    % sum_{h=1}^{k_{l|j}} {n_{h|j,l}/(n_{l|j}+alpha_omega) F(x|omega^{*})}
                    sum_h_mc = sum_h_mc + (nhjl_mc * (nlj_mc+alpha_omega)) .* exp(logprodx1_mc +  logprodx2_mc);
                end % for h
                
                
                % F(t|x,\psi^{*})
                
                tmpcoeff = reshape(gamma_coeff_mat(jj,l,:),size(gamma_coeff0));
                
                prob_t = expit(X_matrix_mc * tmpcoeff');
%                 tmpht1x = binopdf(ones(size(X_mc,1),1),1,prob_t);
                tmpht1x = prob_t;
                wwt1_mc(:,l) = (nlj_mc/(nj_mc + alpha_psi)) .* tmpht1x;
                wwt1_mc(:,l) = wwt1_mc(:,l) .* ( (alpha_omega/(alpha_omega+nlj_mc)) .* f0Xt1_mc + sum_h_mc);
%                 tmpht0x = binopdf(zeros(size(X_mc,1),1),1,prob_t);
                tmpht0x = 1 - prob_t;
                wwt0_mc(:,l) = (nlj_mc/(nj_mc + alpha_psi)) .* tmpht0x;
                wwt0_mc(:,l) = wwt0_mc(:,l) .* ( (alpha_omega/(alpha_omega+nlj_mc)) .* f0Xt0_mc + sum_h_mc);
                sumt1_l_mc = sumt1_l_mc + wwt1_mc(:,l);
                sumt0_l_mc = sumt0_l_mc + wwt0_mc(:,l);
            end % for l
            wt1_mc(:,jj) = nj_mc/(alpha_theta+M) * ( alpha_psi/(nj_mc + alpha_psi) .* f0t1_mc .* f0Xt1_mc + sumt1_l_mc);
            wt0_mc(:,jj) = nj_mc/(alpha_theta+M) * ( alpha_psi/(nj_mc + alpha_psi) .* f0t0_mc .* f0Xt0_mc + sumt0_l_mc);
            partt1_mc = partt1_mc + wt1_mc(:,jj) .* (TXt1_matrix_mc * beta_coeff_mc(jj,:)');
            partt0_mc = partt0_mc + wt0_mc(:,jj) .* (TXt0_matrix_mc * beta_coeff_mc(jj,:)');
        end % for jj
        
        % Cauculate sum_{j=1}^{k} w_{j}(t,x)
        % (This is part of the numerator of E(Y|T=t,X=x,theta^{*},psi^{*},omega^{*},z)
        denompart1_mc = sum(wt1_mc,2);
        denompart0_mc = sum(wt0_mc,2);
        
        % Numerator of E(Y|T=t,X=x,theta^{*},psi^{*},omega^{*},z)
        numt1_mc = partt1_mc + (alpha_theta/(alpha_theta+M)) .* f0t1_mc .* f0Xt1_mc .* e0yt1_mc;
        numt0_mc = partt0_mc + (alpha_theta/(alpha_theta+M)) .* f0t0_mc .* f0Xt0_mc .* e0yt0_mc;
        
        % Denominator of E(Y|T=t,X=x,theta^{*},psi^{*},omega^{*},z)
        denomt1_mc = denompart1_mc + (alpha_theta/(alpha_theta+M)) .* f0t1_mc .* f0Xt1_mc;
        denomt0_mc = denompart0_mc + (alpha_theta/(alpha_theta+M)) .* f0t0_mc .* f0Xt0_mc;
        
       new_psi_rd = mean(numt1_mc./denomt1_mc) - mean(numt0_mc./denomt0_mc);            
