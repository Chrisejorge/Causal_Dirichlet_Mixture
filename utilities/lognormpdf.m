function logval = lognormpdf(x, mu, sigma)
% log normal pdf
% Log of the univariate normal distribution with mean mu and
% standard deviation sigma.  The mean and sigma could be scalars or vectors 
% with the sample number of elements as x. x, mu, sigma should be all
% column vectors or row vectors.
% Returns the log of the pdf.

if nargin < 3
  sigma = 1;
end

if nargin < 2
  mu = 0;
end

if nargin < 1
  error('Requires at least one input argument.');
end

logval =  -log(sigma) -log(2*pi)/2 -0.5 * ((x-mu)./sigma).^2;

% end