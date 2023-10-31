function alpha_omega_new = update_alpha_omega(alpha_shape,alpha_rate,z,alpha_current)
% Define function to update alpha_psi using auxliary variables

k = length(unique(z(:,1)));   % num of y clusters

num_of_tclusters = 0;
for jj = 1: k
    num_of_tclusters =  num_of_tclusters + length(unique(z(z(:,1) ==jj,2)));
end

num_of_xclusters =  size(unique(z,'rows'),1);

nljs = zeros(num_of_tclusters,1);
count = 1;
for jj = 1:k
    kj = length(unique(z(z(:,1) ==jj,2)));
    for l = 1:kj
        tmpidx = z(:,1) ==jj & z(:,2) == l;
        nljs(count) = sum(tmpidx);
        count = count + 1;
    end % for l
end % for jj

xx = betarnd( (alpha_current+1)*ones(size(nljs)),nljs);

zz = rand(size(nljs)).*(alpha_current + nljs) < nljs;

gammaa = alpha_shape + num_of_xclusters - sum(zz);
gammab = alpha_rate - sum(log(xx));

alpha_omega_new = gamrnd(gammaa,1/gammab);

% end