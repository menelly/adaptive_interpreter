"""
setup.py for DNModeling

To install:
- For users: pip install .
- For developers: pip install -e .
"""
from setuptools import setup, find_packages

setup(
    name="dnmodeling",
    version="1.0.0",
    author="Ace, Nova, Ren, and Lumen",
    author_email="noreply@noreply.com",
    description="A revolutionary, domain-aware genetics analysis system.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/menelly/sharedworkspace/tree/main/DNModeling",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=[
        'pyBigWig',
        'requests',
        'pandas',
        'numpy',
        'scikit-learn',
        'joblib',
        'matplotlib',
        'seaborn',
        'xgboost',
        'ensembl_rest',
    ],
    entry_points={
        'console_scripts': [
            'dnmodel-trace=DNModeling.trace_variant:main',
        ],
    },
    include_package_data=True,
)
