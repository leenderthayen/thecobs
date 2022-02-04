import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BSGpy",
    version="0.0.1",
    author="Leendert Hayen",
    author_email="lmhayen@ncsu.edu",
    description="Python port of beta spectrum generator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/leenderthayen/bsgpy",
    project_urls={
        "Bug Tracker": "https://github.com/leenderthayen/bsgpy/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "bsgpy"},
    packages=setuptools.find_packages(where="bsgpy"),
    python_requires=">=3.6",
)
