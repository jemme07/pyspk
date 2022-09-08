import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    print(long_description)

setuptools.setup(
    name="pyspk",
    version=0.1,
    description="Python implementation of the SPK model.",
    url="https://github.com/jemme07/SPK",
    author="Jaime Salcido",
    author_email="j.salcidonegrete@ljmu.ac.uk",
    packages=['pyspk'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy", "scipy"],
)