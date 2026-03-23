from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="plantm6a",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A comprehensive toolkit for plant m6A RNA modification analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/plantm6a",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pyyaml>=5.4",
        "pandas>=1.3.0",
    ],
    entry_points={
        "console_scripts": [
            "plantm6a=plantm6a.cli:main",
        ],
    },
)
