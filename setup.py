from setuptools import setup
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='aqme',
  packages=['aqme'],
  package_data={'aqme': ['templates/*']},
  version='1.0',
  license='MIT',
  description='Automated Quantum Mechanical Environments',
  long_description='Automated Quantum Mechanical Environments',
  long_description_content_type='text/markdown',
  author='Shree Sowndarya S. V., Juan V. Alegre Requena',
  author_email='svss@colostate.edu, juanvi89@hotmail.com',
  keywords=['workflows','computational chemistry','conformational sampling','cheminformatics','quantum mechanics','DFT','automated'],
  url = 'https://github.com/jvalegre/aqme',
  download_url = 'https://github.com/jvalegre/aqme/archive/refs/tags/v_1.0.tar.gz',
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
  ],
  install_requires=[],
  python_requires='>=3.0',
  include_package_data=True,
)
