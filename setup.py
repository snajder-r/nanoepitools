import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanoepitools", # Replace with your own username
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
    python_requires='>=3.6',
)
