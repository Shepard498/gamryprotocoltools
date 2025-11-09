function v = str2num_locale(tokens)
  % Convert a cellstr row of tokens into doubles, flipping comma decimals.
  if ischar(tokens), tokens = regexp(tokens, '[\t,; \s]+', 'split'); end
  v = nan(1, numel(tokens));
  for i = 1:numel(tokens)
    t = tokens{i};
    % Replace decimal comma ONLY when it's between digits
    t = regexprep(t, '(?<=[0-9]),(?=[0-9])', '.');
    v(i) = str2double(t);
  end
end

