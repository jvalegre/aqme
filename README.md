![](Logos/AQME_logo.jpg)
#
## <p align="center"> AQME (Automated Quantum Mechanical Environments)</p>



[![CircleCI](https://img.shields.io/circleci/build/github/jvalegre/aqme?label=Circle%20CI&logo=circleci)](https://app.circleci.com/pipelines/github/jvalegre/aqme)
[![Codecov](https://img.shields.io/codecov/c/github/jvalegre/aqme?label=Codecov&logo=codecov)](https://codecov.io/gh/jvalegre/aqme)
[![Downloads](https://img.shields.io/pepy/dt/aqme?label=Downloads&logo=pypi)](https://www.pepy.tech/projects/aqme)
[![Read the Docs](https://img.shields.io/readthedocs/aqme?label=Read%20the%20Docs&logo=readthedocs)](https://aqme.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/aqme)](https://pypi.org/project/aqme/)

## Documentation  
Full documentation with installation instructions, technical details and examples can be found in [Read the Docs](https://aqme.readthedocs.io).  

Don't miss out the latest hands-on tutorials from our [YouTube channel](https://www.youtube.com/channel/UCHRqI8N61bYxWV9BjbUI4Xw)!  

## Recommended installation
1. (Only once) Create new conda environment: `conda create -n aqme python=3.10`  
2. Activate conda environment: `conda activate aqme`  
3. Install AQME using pip: `pip install aqme`  
4. Install Open Babel: `conda install -y -c conda-forge openbabel=3.1.1`  

* Inexperienced users should visit the *Users with no Python experience* section in our [Read the Docs](https://aqme.readthedocs.io).

## Update the program
1. Update to the latest version: `pip install aqme --upgrade`  

## Developers and help desk
List of main developers and contact emails:  
  - [ ] [Juan V. Alegre-Requena](https://orcid.org/0000-0002-0769-7168), main developer of the CSEARCH, QCORR, QPREP and QDESCP modules. Contact: [jv.alegre@csic.es](mailto:jv.alegre@csic.es)  
  - [ ] [Shree Sowndarya S. V.](https://orcid.org/0000-0002-4568-5854), main developer of the CSEARCH, CMIN and QDESCP modules. Contact: [svss@colostate.edu](mailto:svss@colostate.edu)  
  - [ ] [Brenda Manzanilla](https://orcid.org/0000-0001-5955-6079), developer of the QDESCP module. Contact: [iqmanzanilla@gmail.com](mailto:iqmanzanilla@gmail.com)  
  - [ ] [Turki Alturaifi](https://www.chem.pitt.edu/person/turki-alturaifi), worked in benchmarking the parameters for RDKit-based conformer generation. Contact: [tma53@pitt.edu](mailto:tma53@pitt.edu)  
  - [ ] [Raúl Pérez-Soto](https://orcid.org/0000-0002-6237-2155), worked in refactoring the code and creating the documentation. Contact: [Raul.Perez_Soto@colostate.edu](mailto:Raul.Perez_Soto@colostate.edu)  
  - [ ] [Robert S. Paton](https://orcid.org/0000-0002-0104-4166), research group supervisor and code advisor. Contact: [robert.paton@colostate.edu](mailto:robert.paton@colostate.edu)  

For suggestions and improvements of the code (greatly appreciated!), please reach out through the issues and pull requests options of Github.  

## License
AQME is freely available under an [MIT](https://opensource.org/licenses/MIT) License  

## Reference
If you use any of the AQME modules, please include this citation:  
Alegre-Requena, J. V.; Sowndarya, S.; Pérez-Soto, R.; Alturaifi, T.; Paton, R. AQME: Automated Quantum Mechanical Environments for Researchers and Educators. *Wiley Interdiscip. Rev. Comput. Mol. Sci.* **2023**, *13*, e1663. (DOI: 10.1002/wcms.1663).  
  
Additionally, please include the corresponding references for the following programs:  
  * If you used CSEARCH with RDKit methods or from SMILES: [RDKit](https://www.rdkit.org)  
  * If you used CSEARCH with CREST methods: [CREST](https://crest-lab.github.io/crest-docs)  
  * If you used CMIN with xTB: [xTB](https://xtb-docs.readthedocs.io/en/latest/contents.html)  
  * If you used CMIN with ANI: [ANI](https://github.com/isayev/ASE_ANI)  
  * If you used QCORR: [cclib](https://cclib.github.io/)  
  * If you used QDESCP with xTB: [xTB](https://xtb-docs.readthedocs.io/en/latest/contents.html)
