function logval = logbernpdf(x,pi)
% log binomial pdf
% Log of the  binomial distribution. x and pi should all be row-based matrix or
% columan-based matrix.
% Returns the log of the pdf.

if nargin < 2
  error('Requires at least one input argument.');
end

logval = log(x.*pi + (1-x).*(1-pi) + eps);

% end