function S = concat_cycle(file_list, importer_opts, importer_handle, opts)
%CONCAT_CYCLE  Stitch files into a single cycle-local trace (t starts at 0).
%   S fields: t_cycle, t_local, V, I, fname (cellstr), idx_file (numeric)
%   importer_handle: function handle like @(fp,io) import_gamry_dta(fp, io)

  if nargin < 4, opts = struct(); end
  if nargin < 3 || isempty(importer_handle)
    importer_handle = @import_gamry_dta;
  end
  t_cycle = []; t_local = []; V_all = []; I_all = [];
  fname = {}; idx_file = [];
  t_offset = 0;

  for ii = 1:numel(file_list)
    [t,V,I] = importer_handle(file_list(ii).full, importer_opts);
    tl = t(:) - t(1);
    if ii == 1
      tc = tl;
    else
      tc = t_offset + tl;
    end
    t_offset = tc(end);

    n = numel(tl);
    t_cycle = [t_cycle; tc]; %#ok<AGROW>
    t_local = [t_local; tl]; %#ok<AGROW>
    V_all   = [V_all; V(:)]; %#ok<AGROW>
    I_all   = [I_all; I(:)]; %#ok<AGROW>
    fname   = [fname; repmat({file_list(ii).name}, n, 1)]; %#ok<AGROW>
    idx_file= [idx_file; ii*ones(n,1)]; %#ok<AGROW>
  end

  S = struct('t_cycle', t_cycle, 't_local', t_local, ...
             'V', V_all, 'I', I_all, 'fname', {fname}, 'idx_file', idx_file);
end

