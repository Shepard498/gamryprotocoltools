function progress(msg, frac, opts)
%PROGRESS Lightweight progress logger with optional callback.
%   progress('message', 0.42, opts)
%   If opts.progress_cb exists, calls it; else prints to console.

  try
    if nargin < 2 || isempty(frac), frac = []; end
    if nargin < 3, opts = struct(); end

    if isstruct(opts) && isfield(opts,'progress_cb') && ~isempty(opts.progress_cb)
      opts.progress_cb(msg, frac);
      return;
    end

    if isempty(frac)
      fprintf('[%s] %s\n', caller_tag(), msg);
    else
      fprintf('[%s] %s  (%.0f%%)\n', caller_tag(), msg, 100*max(0,min(1,frac)));
    end
  catch
    % be silent on any logging error
  end
end

function tag = caller_tag()
  % infer a short tag from the caller name, default 'PIPE'
  st = dbstack(2, '-completenames');
  if isempty(st), tag = 'PIPE'; return; end
  nm = st(1).name;
  % collapse long names
  tag = upper(regexprep(nm,'^.*?pipeline_','PIPE_'));
end
