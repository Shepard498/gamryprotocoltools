function R = pipeline_eis(folder, opts)
% PIPELINE_EIS  Stub entry-point for EIS processing.
%   R = pipeline_eis(folder, opts)
% Expected:
%   - Read *.DTA from data/raw (or 'folder' if given)
%   - Export clean tables into data/processed
%   - Save plots into plots/
%
% Notes:
%   - Keep library utilities in octave/lib
%   - Add more pipelines in octave/pipelines (e.g., pipeline_activacion, pipeline_polarizacion)
%
if nargin < 1 || isempty(folder), folder = fullfile(pwd, 'data', 'raw'); end
if nargin < 2, opts = struct(); end

% Example: just list files for now
files = dir(fullfile(folder, '*.DTA'));
fprintf('Found %d DTA files in %s\n', numel(files), folder);

% Return a small report struct
R = struct('n_files', numel(files), 'folder', folder);
end
