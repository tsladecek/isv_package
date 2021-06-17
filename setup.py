"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

CURRENT_DIR = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (CURRENT_DIR / 'README.md').read_text(encoding='utf-8')

# Setting up
setup(
    name="isv",
    version="VERSION",
    author="Tomas Sladecek",
    author_email="tomas.sladecek@geneton.sk",
    description='Automated Interpretation of Structural Copy Number Variants',
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.6, <4',
    install_requires=["numpy", "xgboost", "pandas", "shap", "sklearn-json"],
    keywords=['python', 'machine learning', 'copy number variation'],
    license_files=("LICENSE.txt"),
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        
        'License :: Free for non-commercial use',

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
