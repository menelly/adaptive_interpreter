"""
setup.py for AdaptiveInterpreter

To install:
- For users: pip install .
- For developers: pip install -e .
"""
from setuptools import setup, find_packages

setup(
    name="adaptiveinterpreter",
    version="1.0.0",
    author="Ace, Nova, Ren, and Lumen",
    author_email="shalia@chaoscodex.app",
    description="A revolutionary, domain-aware genetics analysis system.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/menelly/adaptive_interpreter",
    download_url="https://github.com/menelly/adaptive_interpreter/archive/refs/tags/v1.0.0.tar.gz",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
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
            'adaptive-trace=AdaptiveInterpreter.trace_variant:main',
        ],
    },
    include_package_data=True,
)
