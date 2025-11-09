function gui_run()
  % GUI launcher for Gamry pipelines
  % - Pipelines: activacion (ready), polarizacion, ocp, eis (placeholders)
  % - Single input folder
  % - Output file in same folder (auto name; editable)
  % - Alpha textbox (0..1) used only by 'activacion'
  % - Export Excel toggle
  % - Open output folder button after success
  % - Remembers last folder in a small MAT file

  addpath(fullfile(pwd,'lib'));
  addpath(fullfile(pwd,'pipelines'));

  prefs_file = fullfile(pwd, '.gamry_gui_prefs.mat');
  last_folder = pwd;
  if exist(prefs_file, 'file')
    try
      S = load(prefs_file);
      if isfield(S,'last_folder') && ischar(S.last_folder) && exist(S.last_folder,'dir')
        last_folder = S.last_folder;
      end
    catch
      % ignore reading errors
    end
  end

  pipelines = {'activacion', 'polarizacion', 'ocp', 'eis'};

  % ----- figure -----
  f = figure('Name','Gamry Pipelines', 'MenuBar','none', 'NumberTitle','off', ...
             'Position',[200 200 580 230], 'Resize','off');

  % Folder row
  uicontrol(f, 'Style','text', 'String','Folder:', 'HorizontalAlignment','left', ...
               'Position',[20 185 60 20]);
  hFolder = uicontrol(f, 'Style','edit', 'String', last_folder, ...
               'Position',[80 182 390 25]);
  uicontrol(f, 'Style','pushbutton', 'String','Browse...', ...
               'Position',[480 182 80 25], 'Callback', @on_browse);

  % Pipeline row
  uicontrol(f, 'Style','text', 'String','Pipeline:', 'HorizontalAlignment','left', ...
               'Position',[20 145 60 20]);
  hPipe = uicontrol(f, 'Style','popupmenu', 'String', pipelines, ...
               'Position',[80 142 160 25], 'Callback', @on_pipe_change);

  % Alpha row (only for "activacion")
  uicontrol(f, 'Style','text', 'String','Alpha (0..1):', 'HorizontalAlignment','left', ...
               'Position',[260 145 90 20]);
  hAlpha = uicontrol(f, 'Style','edit', 'String','1.0', ...
               'Position',[350 142 60 25]);

  % Export toggle
  hExport = uicontrol(f, 'Style','checkbox', 'String','Export Excel', 'Value',1, ...
               'Position',[430 143 120 22]);

  % Output file row
  uicontrol(f, 'Style','text', 'String','Output name:', 'HorizontalAlignment','left', ...
               'Position',[20 105 90 20]);
  default_name = suggest_outname(pipelines{1});
  hOut = uicontrol(f, 'Style','edit', 'String', default_name, ...
               'Position',[110 102 450 25]);

  % Buttons row
  hRun  = uicontrol(f, 'Style','pushbutton', 'String','Run', 'FontWeight','bold', ...
               'Position',[170 45 110 35], 'Callback', @on_run);
  hOpen = uicontrol(f, 'Style','pushbutton', 'String','Open output folder', ...
               'Position',[300 45 150 35], 'Callback', @on_open_folder, 'Enable','off');

  % Initialize alpha enabled state
  set_alpha_enabled();

  % ------- Callbacks -------
  function on_browse(~,~)
    d = uigetdir(get(hFolder,'String'), 'Seleccione la carpeta con .DTA');
    if ischar(d) && d ~= 0
      set(hFolder,'String', d);
      % update default output name based on pipeline
      p = get_current_pipeline();
      set(hOut, 'String', fullfile(d, suggest_outname(p)));
      % persist
      last_folder = d;
      try, save(prefs_file, 'last_folder'); end
    end
  end

  function on_pipe_change(~,~)
    % enable/disable alpha if needed
    set_alpha_enabled();
    % refresh default output name for current folder
    folder = get(hFolder,'String');
    p = get_current_pipeline();
    set(hOut, 'String', fullfile(folder, suggest_outname(p)));
  end

  function set_alpha_enabled()
    p = get_current_pipeline();
    if strcmpi(p,'activacion')||strcmpi(p,'polarizacion')
      set(hAlpha,'Enable','on');
    else
      set(hAlpha,'Enable','off');
    end
  end

  function p = get_current_pipeline()
    val = get(hPipe,'Value');
    items = get(hPipe,'String');
    if ischar(items), items = cellstr(items); end
    p = items{val};
  end

  function on_run(~,~)
    folder = get(hFolder,'String');
    if ~exist(folder,'dir')
      errordlg('Invalid folder.','Error'); return;
    end
    p = get_current_pipeline();

    % alpha (only for activacion)
    alpha = str2double(get(hAlpha,'String'));
    if isnan(alpha), alpha = 1.0; end
    alpha = max(0, min(1, alpha));

    % export toggle
    do_export = get(hExport,'Value') ~= 0;

    % output filename
    outname = get(hOut,'String');
    if isempty(outname)
      outname = fullfile(folder, suggest_outname(p));
      set(hOut,'String', outname);
    end

    % build opts for pipeline
    opts = struct();
    opts.do_export       = do_export;
    opts.export_filename = outname;

    if strcmpi(p,'activacion')||strcmpi(p,'polarizacion')
      opts.alpha = alpha;
    end

    % Save last folder preference
    last_folder = folder;
    try, save(prefs_file, 'last_folder'); end

    % Run it
    try
      % All actual console logging remains visible in the terminal
      run_pipeline(p, folder, opts);
      % enable "Open output folder"
      if exist(outname,'file') || exist(replace_ext(outname,'.ods'),'file')
        set(hOpen,'Enable','on');
      else
        % even if no file (plots-only), still enable: open folder
        set(hOpen,'Enable','on');
      end
      msgbox('Done.','OK','modal');
    catch err
      errordlg(sprintf('Failed: %s', err.message),'Error');
    end
  end

  function on_open_folder(~,~)
    folder = get(hFolder,'String');
    if ispc()
      system(sprintf('explorer "%s"', folder));
    elseif ismac()
      system(sprintf('open "%s"', folder));
    else
      system(sprintf('xdg-open "%s"', folder));
    end
  end
end

% --------- helpers (file-local) ----------
function s = suggest_outname(pipeline)
  ts = datestr(now,'yyyymmdd_HHMM');
  s = sprintf('Gamry_%s_%s.xlsx', lower(pipeline), ts);
end

function y = replace_ext(x, newext)
  [a,b,~] = fileparts(x);
  y = fullfile(a, [b newext]);
end

