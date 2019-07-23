import setuptools
from codecs import open # To open the README file with proper encoding


# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename, encoding='utf-8') as f:
        return f.read()

setuptools.setup(
    name="arithmat",
    version=readfile("VERSION").strip(),
    author='Giovanni Paolini, Roberto Pagaria',
    author_email='giovanni.paolini@sns.it', # choose a main contact email
    description='Sage implementation of arithmetic matroids and toric arrangements',
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url='https://github.com/giove91/arithmat',
    license='GPLv3', # This should be consistent with the LICENCE file
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Topic :: Software Development :: Build Tools',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 2.7',
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords="SageMath packaging",
    project_urls={
        # 'Documentation': 'https://packaging.python.org/tutorials/distributing-packages/',
        'Source': 'https://github.com/giove91/arithmat',
        'Tracker': 'https://github.com/giove91/arithmat/issues',
    },
    packages=['arithmat'],
    setup_requires=['sage-package'],
    install_requires=['sage-package', 'networkx~=2.2'],
    python_requires='>=2.7, <3',
)
