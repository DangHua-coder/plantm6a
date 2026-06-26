from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="plantm6a",
    version="0.1.0",
    author="PlantM6A contributors",
    description="A toolkit for plant m6A RNA modification analysis",
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
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "pyyaml>=5.4",
    ],
    extras_require={
        "motif": ["pysam>=0.19.0"],
        "conservation": ["parasail>=1.3.0"],
        "visualization": [
            "matplotlib>=3.5.0",
            "scipy>=1.7.0",
            "seaborn>=0.11.0",
            "scikit-learn>=1.0.0",
        ],
        "dev": ["pytest>=7.0.0"],
    },
    entry_points={
        "console_scripts": [
            "plantm6a=plantm6a.cli:main",
        ],
    },
)
