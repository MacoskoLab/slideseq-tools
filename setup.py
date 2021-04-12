#!/usr/bin/env python

import io
import glob
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="slideseq-tools",
    version="0.1",
    license="MIT License",
    description="Scripts processing slideseq data",
    long_description=read("README.md"),
    url="https://github.com/MacoskoLab/slideseq-tools",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("src/*.py")
    ],
    include_package_data=True,
    zip_safe=False,
    install_requires=["Click", "numpy", "pandas", "plotnine", "matplotlib"],
    entry_points={"console_scripts": ["submit_job = slideseq.submit_job:main"]},
)
