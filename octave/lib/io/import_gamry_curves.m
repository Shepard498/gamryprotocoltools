function R = import_gamry_curves(filepath)
  T = parse_dta_table(filepath, 'CURVE');
  fn = fieldnames(T);
  get = @(names) pick_field(T, fn, names);
  t = get({'time_s','time','t_s','t'});
  v = get({'potential_v','ewe_v','v_v','v','vdiff_v'});
  i = get({'current_a','i_a','i'});
  if isempty(t) || isempty(v) || isempty(i)
    error('Unrecognized CURVE columns in %s', filepath);
  end
  R = struct('t_s', t(:), 'V_V', v(:), 'I_A', i(:));
end

