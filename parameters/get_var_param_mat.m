function x_var_param_mat = get_var_param_mat(z,k,kjs,kljs,x_var_param)
% Get matrix for x_var

maxNumTclusters = max(z(:,2));
maxNumXclusters = max(z(:, 3));

x_var_param_mat = zeros(k+1,maxNumTclusters+1,maxNumXclusters+1,size(x_var_param,2));

count = 1;
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        klj = kljs(jj,l);
        for h = 1:klj
            x_var_param_mat(jj,l,h,:) = x_var_param(count,:);
            count = count + 1;
        end % for h
    end % for l
end % for jj

% end
