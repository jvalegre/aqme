from setuptools import setup
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='DBGEN',
  packages=['DBGEN'],
  version='1.0',
  description='Conformer generation followed DFT analysis',
  long_description=long_description,
  long_description_content_type='text/markdown',
  author='',
  author_email='robert.paton@colostate.edu',
  keywords=['compchem', 'informatics'],
  classifiers=[],
  install_requires=["numpy","pandas"],
  python_requires='>=3.0',
  include_package_data=True,
)
