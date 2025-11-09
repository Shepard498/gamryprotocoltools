function export_excel(xlsx_path, sheets)
  pkg('load', 'io');
  for i = 1:size(sheets,1)
    sname = sanitize_sheetname(sheets{i,1});
    cells = sheets{i,2};
    xlswrite(xlsx_path, cells, sname);
  end
end
