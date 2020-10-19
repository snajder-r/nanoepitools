import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanoepitools",
    version="0.1.0",
    author="Rene Snajder",
    author_email="r.snajder@dkfz-heidelberg.de",
    description="Package for epigenetic analyses from Nanopore data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/snajder-r/nanoepitools´",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        "numpy>=1.14.0",
        "h5py>=2.9.0",
        "pandas>=1.0.3",
        "matplotlib",
        ],
    python_requires='>=3.6',
)
