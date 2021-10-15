from setuptools import setup, find_packages

path_to_readme = None

if path_to_readme is not None:
    with open(path_to_readme, "r") as desc_file:
        long_description = desc_file.read()
else:
    long_description = ''

path_to_version_file = "./structman/_version.py"

with open(path_to_version_file) as version_file:
    exec(version_file.read().strip())

setup(
    name="StructMAn",
    version=__version__,
    description="Structural Mutation Annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='LGPL-2.1',
    author="Alexander Gress",
    maintainer="Alexander Gress",

    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    include_package_data=True,
    install_requires=[
        "biopython>=1.78",
        "python-igraph>=0.8.3",
        "matplotlib>=3.3.2",
        "numpy>=1.19.2",
        "psutil>=5.8.0",
        "pymysql>=1.0.2",
        "python-igraph>=0.8.3",
        "ray==1.6.0",
        "msgpack==1.0.2",
        "zstd==1.5.0.2"
    ],
    python_requires=">=3.8, <4",
    keywords="bioinformatics",
    entry_points={
        "console_scripts": ["structman = structman.structman_main:structman_cli"],
    },
)
