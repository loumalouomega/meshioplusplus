function mio_bar_chart(labels, values, chart_title, color)
%MIO_BAR_CHART  Horizontal bar chart: file sizes per format, cell counts per type.
%   mio_bar_chart(LABELS, VALUES, TITLE) draws one bar per label, the first
%   label at the top, in a new figure. LABELS is a cellstr. Notebook helper,
%   not part of the package API.

    if nargin < 4
        color = [59 130 246] / 255;
    end
    labels = labels(:);
    values = double(values(:));
    n = numel(values);
    figure();
    barh(1:n, flipud(values), 'facecolor', color, 'edgecolor', 'none');
    set(gca, 'ytick', 1:n, 'yticklabel', flipud(labels), 'ylim', [0.4, n + 0.6]);
    title(chart_title, 'interpreter', 'none');
    grid on
end
