function alpha_psi_new = update_alpha_psi(alpha_shape,alpha_rate,z,alpha_current)
% Define function to update alpha_psi using auxliary variables

k = length(unique(z(:,1)));

num_of_tclusters = 0;
for jj = 1: k
    num_of_tclusters =  num_of_tclusters + length(unique(z(z(:,1) == jj,2)));
end

% njs = zeros(k,1);
% for jj = 1:k
%     njs(jj) = sum(z(:,1) == jj);
% end

njs = crosstab(z(:,1));

xx = betarnd( (alpha_current+1)*ones(size(njs)),njs);

zz = rand(size(njs)).*(alpha_current + njs) < njs;

gammaa = alpha_shape + num_of_tclusters - sum(zz);
gammab = alpha_rate - sum(log(xx));

alpha_psi_new = gamrnd(gammaa,1/gammab);

% end