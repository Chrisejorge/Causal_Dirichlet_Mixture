function y = logdet(A)
% log det A

U = chol(A);
y = 2*sum(log(diag(U)));

% end