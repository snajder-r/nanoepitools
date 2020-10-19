from pathlib import Path

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bokeh.plotting import figure, output_file


def set_if_not_in_dict(d, k, v):
    if k not in d.keys():
        d[k] = v


def set_figure_defaults(fig):
    fig.set_tight_layout(True)
    fig.patch.set_facecolor("w")
    fig.autolayout = False


def density_plot(x, y, *argc, **kwargs):
    density = scipy.stats.gaussian_kde(y)
    plt.plot(density(x), *argc, **kwargs)


def plot_2d_density(x, y, nbins=50, cmap=plt.cm.BuGn_r):
    k = scipy.stats.gaussian_kde((x, y))
    xi, yi = np.mgrid[x.min() : x.max() : nbins * 1j, y.min() : y.max() : nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading="gouraud", cmap=cmap)
    plt.contour(xi, yi, zi.reshape(xi.shape))


class PDFPagesWrapper:
    def __init__(self, pa, *argv):
        self.pa = pa
        self.pdf = PdfPages(*argv)

    def __enter__(self):
        self.pdf.__enter__()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pa.pdf = None
        self.pdf.__exit__(exc_type, exc_val, exc_tb)


class PlotArchiver:
    def __init__(self, project, config=None):
        self.project = project
        if config is None:
            config = dict()
        set_if_not_in_dict(
            config, "plot_archive_dir", Path.home().joinpath("nanoepitools_plots")
        )
        set_if_not_in_dict(config, "filetype", "pdf")
        self.config = config

        self.project_path = Path(self.config["plot_archive_dir"]).joinpath(self.project)
        self.pdf: PDFPagesWrapper = None

    def ensure_project_path_exists(self):
        self.project_path.mkdir(parents=True, exist_ok=True)

    def get_plot_path(self, key, filetype=None):
        if filetype is None:
            filetype = self.config["filetype"]
        filename = "{key}.{ft}".format(key=key, ft=filetype)
        return self.project_path.joinpath(filename)

    def savefig(self, key="figure", fig=None, close=True):
        self.ensure_project_path_exists()
        if fig is None:
            fig = plt.gcf()
        if self.pdf is not None:
            self.pdf.pdf.savefig(fig)
        else:
            path = self.get_plot_path(key)
            fig.savefig(path)
        if close:
            plt.close(fig)

    def saveandshow(self, key="figure", fig=None):
        if fig is None:
            fig = plt.gcf()
        self.savefig(key, fig=fig, close=False)
        fig.show()
        # This is a workaround needed to ensure the figure is actually
        # displayed in jupyter notebook. If we don't wait for a bit and close
        # it right away, it may not be displayed
        plt.pause(0.0001)
        plt.close(fig)

    def bokeh_open_html(self, key):
        self.ensure_project_path_exists()
        path = self.get_plot_path(key, "html")
        output_file(path)

    def open_multipage_pdf(self, key):
        self.ensure_project_path_exists()
        path = self.get_plot_path(key, "pdf")
        self.pdf = PDFPagesWrapper(self, path)
        return self.pdf

    def subplots(self, *args, **kwargs):
        if "dpi" not in kwargs:
            kwargs["dpi"] = 200
        fig, axes = plt.subplots(*args, **kwargs)
        set_figure_defaults(fig)
        return fig, axes

    def figure(self, **kwargs):
        """
        simply calls matplotlib.pyplot.figure() but then sets some default
        parameters that I find handy
        :param argc: parameters for plt.figure
        :returns: return value if plt.figure
        """
        if "dpi" not in kwargs:
            kwargs["dpi"] = 200
        fig = plt.figure(**kwargs)
        set_figure_defaults(fig)
        return fig
