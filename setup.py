"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import os
import pathlib

CURRENT_DIR = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (CURRENT_DIR / 'README.md').read_text(encoding='utf-8')

VERSION = open(os.path.join(CURRENT_DIR, 'VERSION')).read().strip()

# Setting up
setup(
    name="isv",
    version=VERSION,
    author="Tomas Sladecek",
    author_email="tomas.sladecek@geneton.sk",
    description='Automated Interpretation of Structural Copy Number Variants',
    long_description_content_type="text/markdown",
    long_description=long_description,
    package_dir={'': 'isv'},
    packages=find_packages(where='isv'),
    python_requires='>=3.6, <4',
    install_requires=["numpy>=1.20.0",
                      "xgboost>=1.4.0",
                      "pandas>=1.2.0",
                      "shap>=0.39.0",
                      "sklearn-json>=0.1.0"
                      ],
    keywords=['python', 'machine learning', 'copy number variation'],
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
)
