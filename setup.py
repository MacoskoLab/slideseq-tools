#!/usr/bin/env python

import glob
import io
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="slideseq-tools",
    version="0.2.1",
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
    install_requires=[
        "Click",
        "PyYAML",
        "google-api-python-client",
        "google-auth",
        "google-cloud-secret-manager",
        "matplotlib",
        "networkx",
        "numpy",
        "openpyxl",
        "pandas",
        "plotnine",
        "scikit-learn",
        "scipy",
    ],
    extras_require={"dev": ["black", "isort", "flake8", "pre-commit"]},
    entry_points={
        "console_scripts": [
            "submit_flowcell = slideseq.pipeline.submit_slideseq:main",
            "align_library = slideseq.pipeline.alignment:main",
            "process_library = slideseq.pipeline.processing:main",
            "build_ref = slideseq.pipeline.reference:main",
        ]
    },
)
