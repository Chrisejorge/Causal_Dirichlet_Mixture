function x_pi_param_mat = get_pi_param_mat(z,k,kjs,kljs,x_pi_param)
% Get matrix for x_pi

maxNumTclusters = max(z(:,2));
maxNumXclusters = max(z(:, 3));

x_pi_param_mat = zeros(k+1,maxNumTclusters+1,maxNumXclusters+1,size(x_pi_param,2));

count = 1;
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        klj = kljs(jj,l);
        for h = 1:klj
            x_pi_param_mat(jj,l,h,:) = x_pi_param(count,:);
            count = count + 1;
        end % for h
    end % for l
end % for jj
% end