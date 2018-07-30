import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="canopus",
    version="0.0.1",
    author="Kai Dührkop",
    author_email="kai.duehrkop@uni-jena.de",
    description="visualization for CANOPUS with Jupyter notebook",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaibioinfo/canopus_treemap",
    packages=setuptools.find_packages(),
    package_data={'': ['*.js']},
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)