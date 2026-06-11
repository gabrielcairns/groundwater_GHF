% Function to make nice colorbars
function cb = make_nice_colorbar(ax,location,label,varargin)

cb = colorbar(ax,location);

% Add a label, make everything size 14 latex
cb.Label.String = label;
cb.Label.FontSize = 14;
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

% Optional argument to modify limits from default
if nargin >= 4
    lims = varargin{1};
    cb.Limits = lims;
end

end