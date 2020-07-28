import unittest
from pathlib import Path

import matplotlib.pyplot as plt

from nanoepitools.plotting.general_plotting import PlotArchiver


def expected_default_folder():
    return Path.home().joinpath('nanoepitools_plots').joinpath(
        'testproject')


def clean_if_exists(filepath: Path):
    if filepath.exists():
        filepath.unlink()


class TestPlotArchive(unittest.TestCase):

    def test_default_path(self):
        pa = PlotArchiver('testproject')
        self.assertEqual(expected_default_folder().joinpath(
            'my_test_plot.pdf'), pa.get_plot_path('my_test_plot'))

    def test_custom_path(self):
        config = {'plot_archive_dir': '/tmp/test'}
        pa = PlotArchiver('testproject', config=config)
        self.assertEqual(Path('/tmp/test').joinpath('testproject').joinpath(
            'my_test_plot.pdf'), pa.get_plot_path('my_test_plot'))

    def test_create_plot(self):
        expected_filepath = expected_default_folder().joinpath('testplot.pdf')
        clean_if_exists(expected_filepath)

        plt.plot([0, 1], [0, 1])
        pa = PlotArchiver('testproject')
        pa.savefig('testplot', plt.gcf())

        self.addCleanup(expected_filepath.unlink)
        self.assertTrue(expected_filepath.exists(), 'File does not exist')
