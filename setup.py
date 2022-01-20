from setuptools import setup, find_packages
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='aqme',
  packages=find_packages(exclude=['tests']),
  package_data={'aqme': ['templates/*']},
  version='0.0',
  license='MIT',
  description='Automated Quantum Mechanical Environments',
  long_description='Automated Quantum Mechanical Environments',
  long_description_content_type='text/markdown',
  author='Shree Sowndarya S. V., Juan V. Alegre Requena',
  author_email='svss@colostate.edu, juanvi89@hotmail.com',
  keywords=['workflows','computational chemistry','conformational sampling','cheminformatics','quantum mechanics','DFT','automated'],
  url = 'https://github.com/jvalegre/aqme',
  download_url = 'https://github.com/jvalegre/aqme/archive/refs/tags/0.0.tar.gz',
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',      #Specify which python versions you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
  ],
  install_requires=['PyYAML>=5.3',
    'pandas>=1.0.3',
    'progress>=1.5',
    'ase>=3.19.1',
    'numpy>=1.18.4',
    'cclib>=1.6.3',
    'matplotlib>=3.2.2',
    'seaborn>=0.10.1',
    'torch>=1.5.0+cpu',
    'torchvision>=0.6.0+cpu',
    'torchani>=2.2'],
  python_requires='>=3.0',
  include_package_data=True,
)
