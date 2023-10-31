function newval = log1m( x )
% Log of 1-x, more accurate than log(1-x)
% Copyright (C) Adi Lin

newval = log1p(-x);

% end

