function mio_histogram(values, chart_title, xlabel_text, color)
%MIO_HISTOGRAM  Histogram of a numeric vector, in a new figure.
%   mio_histogram(VALUES, TITLE, XLABEL) bins VALUES into 12 bins. Notebook
%   helper, not part of the package API.

    if nargin < 4
        color = [59 130 246] / 255;
    end
    figure();
    [counts, centres] = hist(double(values(:)), 12);
    bar(centres, counts, 1, 'facecolor', color, 'edgecolor', [1 1 1]);
    title(chart_title, 'interpreter', 'none');
    xlabel(xlabel_text, 'interpreter', 'none');
    ylabel('count');
end
