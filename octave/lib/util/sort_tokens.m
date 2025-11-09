function S = sort_tokens(files, mode)
  % Sort according to protocol tokens in filenames.
  % mode = 'activation' or 'polarization'
  % Expected patterns:
  %  Activation: Activacion_(Asc|Dsc)_(proto)_#CYCLE_#FILE.DTA
  %  Polarization: Curva_Polarizacion_(Asc|Dsc)_(proto)_#FILE.DTA
  rec = struct('name', '', 'path', '', 'dir', '', 'proto', '', 'dirflag', 0, 'cycle', NaN, 'file', NaN);
  R = repmat(rec, numel(files), 1);
  for k = 1:numel(files)
    f = files(k);
    R(k).name = f.name; R(k).path = fullfile(f.folder, f.name); R(k).dir = f.folder;
    s = f.name;
    switch mode
      case 'activation'
        m = regexp(s, '^Activaci[oó]n_(Asc|Dsc)_([^_]+)_#(\d+)_#(\d+)\.', 'tokens', 'once', 'ignorecase');
        if ~isempty(m)
          R(k).dirflag = strcmpi(m{1}, 'Asc');
          R(k).proto = m{2};
          R(k).cycle = str2double(m{3});
          R(k).file  = str2double(m{4});
        end
      case 'polarization'
        m = regexp(s, '^Curva_Polarizaci[oó]n_(Asc|Dsc)_([^_]+)_#(\d+)\.', 'tokens', 'once', 'ignorecase');
        if ~isempty(m)
          R(k).dirflag = strcmpi(m{1}, 'Asc');
          R(k).proto = m{2};
          R(k).file  = str2double(m{3});
        end
    end
  end
  switch mode
    case 'activation'
      key = [[R.cycle]' [~[R.dirflag]]' [R.file]'];
      [~, idx] = sortrows(key, [1 2 3]);
    case 'polarization'
      key = [[~[R.dirflag]]' [R.file]'];
      [~, idx] = sortrows(key, [1 2]);
  end
  S = R(idx);
end

