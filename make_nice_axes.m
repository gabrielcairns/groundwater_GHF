% Function to initialise axes with latex axis labels and optional title
function ax = make_nice_axes(ax,xlabeltext,ylabeltext,varargin)

cla(ax);
hold(ax, 'on');

% Readable-fontsize latex axis labels
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';

% Adding text labels
xlabel(ax, xlabeltext, 'FontSize', 14, 'Interpreter', 'Latex')
ylabel(ax, ylabeltext, 'FontSize', 14, 'Interpreter', 'Latex')

% Optional title
if nargin >= 4
    titletext = varargin{1};
    title(ax, titletext, 'FontSize', 16, 'Interpreter', 'Latex')
end
    
end