import matplotlib.figure
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages


def new_ax(pdf_pages: PdfPages, include_fig=False):
    """
    Simple helper to make a new figure and provide the axes for plotting,
    then save to an open PDF when finished.
    """
    fig = matplotlib.figure.Figure(figsize=(8, 8))
    ax: Axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    if include_fig:
        yield fig, ax
    else:
        yield ax

    pdf_pages.savefig(fig)
