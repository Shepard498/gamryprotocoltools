function v = str2num_locale(tokens)
  % Convert cellstr tokens to double, flipping comma decimals to dot.
  v = nan(size(tokens));
  for i = 1:numel(tokens)
    t = tokens{i};
    t = regexprep(t, '(?<=[0-9]),(?=[0-9])', '.');
    v(i) = str2double(t);
  end
end
