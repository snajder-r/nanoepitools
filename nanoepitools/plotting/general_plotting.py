from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bokeh.plotting import figure, output_file


def set_if_not_in_dict(d, k, v):
    if k not in d.keys():
        d[k] = v


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
        set_if_not_in_dict(config, 'plot_archive_dir',
                           Path.home().joinpath('nanoepitools_plots'))
        set_if_not_in_dict(config, 'filetype', 'pdf')
        self.config = config

        self.project_path = Path(self.config['plot_archive_dir']).joinpath(
            self.project)
        self.pdf: PDFPagesWrapper = None

    def ensure_project_path_exists(self):
        self.project_path.mkdir(parents=True, exist_ok=True)

    def get_plot_path(self, key, filetype=None):
        if filetype is None:
            filetype = self.config['filetype']
        filename = '{key}.{ft}'.format(key=key, ft=filetype)
        return self.project_path.joinpath(filename)

    def savefig(self, key='figure', fig=None, close=True):
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

    def saveandshow(self, key, fig=None):
        if fig is None:
            fig = plt.gcf()
        self.savefig(key, fig=fig, close=False)
        fig.show()
        plt.close(fig)

    def bokeh_open_html(self, key):
        self.ensure_project_path_exists()
        path = self.get_plot_path(key, 'html')
        output_file(path)

    def open_multipage_pdf(self, key):
        self.ensure_project_path_exists()
        path = self.get_plot_path(key, 'pdf')
        self.pdf = PDFPagesWrapper(self, path)
        return self.pdf

