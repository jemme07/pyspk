"""Legacy setuptools entrypoint.

This repository is primarily configured via `pyproject.toml` (PEP 517/518).
`setup.py` is kept for compatibility with older workflows.
"""

from __future__ import annotations

from pathlib import Path

from setuptools import find_packages, setup


def _read_readme() -> str:
    readme_path = Path(__file__).with_name("README.md")
    return readme_path.read_text(encoding="utf-8")


setup(
    name="pyspk",
    version="2.0.0",
    description=(
        "Python package to predict the suppression of the total matter power spectrum "
        "due to baryonic physics"
    ),
    url="https://github.com/jemme07/pyspk",
    author="Jaime Salcido",
    author_email="j.salcidonegrete@ljmu.ac.uk",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    long_description=_read_readme(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy>=1.22", "scipy>=1.8", "pydantic>=2"],
    include_package_data=True,
)
