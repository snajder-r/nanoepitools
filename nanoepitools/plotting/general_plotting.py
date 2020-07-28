import matplotlib.pyplot as plt
from pathlib import Path


def set_if_not_in_dict(d, k, v):
    if k not in d.keys():
        d[k] = v


class PlotArchiver:
    def __init__(self, project, config=None):
        self.project = project
        if config is None:
            config = dict()
        set_if_not_in_dict(config, 'plot_archive_dir', Path.home().joinpath(
            'nanoepitools_plots'))
        set_if_not_in_dict(config, 'filetype', 'pdf')
        self.config = config

        self.project_path = Path(self.config['plot_archive_dir']).joinpath(
            self.project)

    def ensure_project_path_exists(self):
        self.project_path.mkdir(parents=True, exist_ok=True)

    def get_plot_path(self, key):
        filename = '{key}.{ft}'.format(key=key, ft=self.config['filetype'])
        return self.project_path.joinpath(filename)

    def savefig(self, key, fig=None):
        self.ensure_project_path_exists()
        if fig is None:
            fig = plt.gcf()
        path = self.get_plot_path(key)
        fig.savefig(path)

    def saveandshow(self, key, fig=None):
        if fig is None:
            fig = plt.gcf()
        self.savefig(key, fig=fig)
        fig.show()


