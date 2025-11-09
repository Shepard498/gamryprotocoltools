function y = smooth_signal(x, window, method)
%SMOOTH_SIGNAL Unified signal smoothing with different methods
%   y = smooth_signal(x, window, method)
%   method: 'median' or 'mean' (default 'median')
%   window: odd window length (will be forced odd if even)
  
  if nargin < 3 || isempty(method), method = 'median'; end
  if nargin < 2 || isempty(window), window = 5; end
  
  % Ensure odd window
  window = max(1, 2*floor((window-1)/2) + 1);
  
  % Early exit for trivial cases
  x = x(:);
  if window == 1 || numel(x) < 3
    y = x;
    return;
  end
  
  % Apply smoothing
  n = numel(x);
  y = zeros(size(x));
  h = (window-1)/2;
  
  for i = 1:n
    i1 = max(1, i-h);
    i2 = min(n, i+h);
    switch lower(method)
      case 'median'
        y(i) = median(x(i1:i2));
      case 'mean'
        y(i) = mean(x(i1:i2));
      otherwise
        error('Unknown smoothing method: %s', method);
    end
  end
end