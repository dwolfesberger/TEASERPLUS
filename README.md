![E.ON EBC RWTH Aachen University](./pictures/pyside_teco.png)

# Teco - An extension to the TEASER+ tool for environmental impacts with input and output functionalities for CityGML and Energy ADE 

<!---[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)-->

The Teco extension is developed by members of the "Institute of Energy Efficiency and Sustainable Building (e3D), RWTH Aachen University" using Python 3.5+.
This extension is based on TEASER+ as well as the "Tool for Energy Analysis and Simulation for Efficient Retrofit (TEASER)" and can be used to import and export CityGML data sets with or without Energy ADE version 1.0. 
The Teco extension aims to help simulation scientists facilitate the determination of environmental impacts (prominently GWP) using basic input and CityGML models.
This GitHub page will be used to further develop the extension.
<!---and make it available under the [MIT License](https://gitlab.e3d.rwth-aachen.de/e3d-software-tools/citybit/citybit/-/blob/master/License/LICENSE).-->

If you have any questions regarding Teco feel free to contact us at: [schildt@e3d.rwth-aachen.de](mailto:schildt@e3d.rwth-aachen.de).

If you have any questions regarding TEASER+ feel free to contact us at: [shamovich@e3d.rwth-aachen.de](mailto:shamovich@e3d.rwth-aachen.de) or [cityatb@e3d.rwth-aachen.de](mailto:cityatb@e3d.rwth-aachen.de).

## Description

Ambitious sustainability and carbon reduction goals in the building sector require the large-scale determination of environmental impacts. Building energy demands in the operational phase are a crucial element to be considered. To fully reflect the building life cycle, emissions from all life phases - construction, operation, dismantling, re-usage - have to be analysed. This goes in accordance with the 11th United Nations Sustainable Goal of making cities and human settlements sustainable. 
At the same time, this analysis is generally data-intensive. With the increasing availability of 3D building models, particularly CityGML datasets, using such data is advantageous to the UBEM community and the life cycle assessment workflow. Since CityGML data is mainly geometric, Teco offers an enrichment procedure for building materials, utilities and respective environmental indicators based on DIN EN 15804. Thus, the Teco extension enhances the existing feature set and abilities of the TEASER+ tool, and allows the integration of 3D city models for urban-scale, ecological life cycle analyses.

## Version

Teco is currently being developed and is based on the TEASER version 0.7.6. 
The current release enables the addition of buildings in the GUI based on relevant input parameters. The next version will comprise the full usability of CityGML data.

<!---## How to use TEASER+-->


### Dependencies

Teco is currently being developed using Python 3.5+ and PySide2. Beside the basic dependencies from the TEASER tool, Teco uses BuildingsPy version 2.1.0 to automate the simulation process. The current developments of Teco use lxml for importing and exporting CityGML datasets with or without the Energy ADE.

### Installation

Teco can be used by cloning or downloading the whole Teco package from the GIT Repository. In future, the user can to run the "main.py" for loading the GUI.

### How to contribute to the development of Teco

You are invited to contribute to the development of Teco. You may report any issues by sending us an email to [schildt@e3d.rwth-aachen.de](mailto:schildt@e3d.rwth-aachen.de).

## How to cite TEASER+
The following articles can be used to cite TEASER+:
+ Urban energy simulations using open CityGML models: A comparative analysis. Malhotra, A., Shamovich, M., Frisch, J., & van Treeck, C.. Energy and Buildings, 255, 111658, January 2022.

### TEASER+ related publications
+ Parametric Study of the Different Level of Detail of CityGML and Energy-ADE Information for Energy Performance Simulations. Malhotra, A., Shamovich, M., Frisch, J., & van Treeck, C.. In Proc. 16th IBPSA Conf (Vol. 16, pp. 3429-3436), September 2019.


## How to cite original TEASER

+ TEASER: an open tool for urban energy modelling of building stocks. Remmen P., Lauster M., Mans M., Fuchs M., Osterhage T., Müller D.. Journal of Building Performance Simulation, February 2017,
[pdf](http://dx.doi.org/10.1080/19401493.2017.1283539),
[bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite_jbps.bib)

### TEASER related publications

+ CityGML Import and Export for Dynamic Building Performance Simulation in Modelica. Remmen P.,
Lauster M., Mans M., Osterhage T., Müller D.. BSO16, p.329-336, September 2016,
[pdf](http://www.ibpsa.org/proceedings/BSO2016/p1047.pdf),
[bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite.bib)

+ Scalable Design-Driven Parameterization of Reduced Order Models Using Archetype Buildings with TEASER.
Lauster M., Mans M., Remmen P., Fuchs M., Müller D.. BauSIM2016, p.535-542, September 2016,
[pdf](https://www.researchgate.net/profile/Moritz_Lauster/publication/310465372_Scalable_Design-Driven_Parameterization_of_Reduced_Order_Models_using_Archetype_Buildings_with_TEASER/links/582ee96908ae004f74be1fb0.pdf?origin=publication_detail&ev=pub_int_prw_xdl&msrp=eEyK6WYemhC8wK7xkMEPRDO4obE4uxBN4-0BdBy1Ldwhy9FhCe1pXfNObJYubvC_aZN0IWDPf9uayBo3u79bsZvg3hzUoLoYRatES2ARH8c.B2cYwSICt0IOa7lD-4oAiEa_3TtrO-7k-1W9chuNQwr_VNMCpZ5ubSb-eY2D77rGUP4S6wS8m6vudUUbMlXbQQ.Cledgd1Q9fPp11nYGpcpKNhSS6bVTqAEXeMZPkiV3HsJxcVWTFj4Hr_jmLZ0MOzDxbDEZObcGiKfmTL_9k_59A)

+ Refinement of Dynamic Non-Residential Building Archetypes Using Measurement Data and Bayesian Calibration
Remmen P., Schäfer J., Müller D.. Building Simulation 2019, September 2019,
[pdf](https://www.researchgate.net/publication/337925776_Refinement_of_Dynamic_Non-Residential_Building_Archetypes_Using_Measurement_Data_and_Bayesian_Calibration)

+ Selecting statistical indices for calibrating building energy models. Vogt, M., Remmen P., Lauster M., Fuchs M. , Müller D.. Building and Environment 144, pages 94-107, October 2018. [bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite_be.bib)

+  A parametric study of TEASER where all functions and  parameters used in TEASER are gathered and explained. The publication can be found [here](https://publications.rwth-aachen.de/record/749801/files/749801.pdf).

<!---## License

CityBIT is released by RWTH Aachen University, E3D - Institute of Energy Efficiency and Sustainable Building, under the [MIT License](https://gitlab.e3d.rwth-aachen.de/e3d-software-tools/citybit/citybit/-/blob/master/License/LICENSE).-->

## Acknowledgements

The developers of Teco would like to thank the Institute of Energy Efficient Building and Indoor Climate (EBC), E.ON Energy Research Center, RWTH Aachen University for their effort in developing TEASER and for making it available open-source. 




<!---


# TEASER -  Tool for Energy Analysis and Simulation for Efficient Retrofit

[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)
[![Coverage Status](https://coveralls.io/repos/github/RWTH-EBC/TEASER/badge.svg)](https://coveralls.io/github/RWTH-EBC/TEASER)
[![Build Status](https://travis-ci.org/RWTH-EBC/TEASER.svg?branch=master)](https://travis-ci.org/RWTH-EBC/TEASER.svg?branch=master)

TEASER (Tool for Energy Analysis and Simulation for Efficient Retrofit) allows
fast generation of archetype buildings with low input requirements and the
export of individual dynamic simulation models for the below-mentioned Modelica
libraries. These libraries all use the framework of [Modelica IBPSA
library](https://github.com/ibpsa/modelica). TEASER is being developed at the
[RWTH Aachen University, E.ON Energy Research Center, Institute for Energy
Efficient Buildings and Indoor
Climate](https://www.ebc.eonerc.rwth-aachen.de/cms/~dmzz/E-ON-ERC-EBC/?lidx=1).

 * [AixLib](https://github.com/RWTH-EBC/AixLib)
 * [Buildings](https://github.com/lbl-srg/modelica-buildings)
 * [BuildingSystems](https://github.com/UdK-VPT/BuildingSystems)
 * [IDEAS](https://github.com/open-ideas/IDEAS).

The full documentation of TEASER including examples and description of modules,
classes and functions can be found at the website:

 * http://rwth-ebc.github.io/TEASER/

This GitHub page will be used to further develop the package and make it
available under the
[MIT License](https://github.com/RWTH-EBC/TEASER/blob/master/License.md).

If you have any questions regarding TEASER feel free to contact us at
[ebc-teaser@eonerc.rwth-aachen.de](mailto:ebc-teaser@eonerc.rwth-aachen.de).

If you want to use TEASER without installation, you can use out TEASER webtool, which
will generate a Modelica model and provide this as download:

 * [http://teaser.eonerc.rwth-aachen.de](http://teaser.eonerc.rwth-aachen.de)

## Description

Energy supply of buildings in urban context currently undergoes significant
changes. The increase of renewable energy sources for electrical and thermal
energy generation will require flexible and secure energy storage and
distribution systems. To reflect and consider these changes in energy systems
and buildings, dynamic simulation is one key element, in particular when it
comes to thermal energy demand on minutely or hourly scale.
Sparse and limited access to detailed building information as well as computing
times are challenges for building simulation on urban scale. In addition,
data acquisition and modeling for Building Performance Simulation (BPS) are
time consuming and error-prone. To enable the use of BPS on urban scale we
present the TEASER tool, an open framework for urban energy modeling of
building stocks. TEASER provides an easy interface for multiple data sources,
data enrichment (where necessary) and export of ready-to-run Modelica simulation
models for all libraries supporting the
[Modelica IBPSA library](https://github.com/ibpsa/modelica).


## Version

TEASER is a ongoing research project, the current version is 0.7.6, which is
still a pre-release.

## How to use TEASER

### Dependencies

TEASER is currently tested against Python 3.6 and 3.7. Older versions of Python may
still work, but are no longer actively supported.
Using a Python distribution is recommended as they already contain (or easily
support installation of) many Python packages (e.g. SciPy, NumPy, pip, PyQT,
etc.) that are used in the TEASER code. Two examples of those distributions are:

1. https://winpython.github.io/ WinPython comes along with a lot of Python
packages (e.g. SciPy, NumPy, pip, PyQT, etc.)..
2. http://conda.pydata.org/miniconda.html Conda is an open source package
management  system and environment management system for installing multiple
versions of software  packages and their dependencies and switching easily
between them.

In addition, TEASER requires some specific Python packages:

1. Mako: template Engine
  install on a python-enabled command line with `pip install -U mako`
2. pandas: popular data analysis library
  install on a python-enabled command line with `pip install -U pandas`
3. pytest: Unit Tests engine
  install on a python-enabled command line with `pip install -U pytest`

### Installation

The best option to install TEASER is to use pip:

`pip install teaser`

If you actively develop TEASER you can clone this repository by using:

 `git clone [SSH-Key/Https]`

and then run:

 `pip install -e [Path/to/your/Teaser/Clone]` which will install the local version of TEASER.


### How to contribute to the development of TEASER
You are invited to contribute to the development of TEASER. You may report any issues by using the [Issues](https://github.com/RWTH-EBC/TEASER/issues) button.
Furthermore, you are welcome to contribute via [Pull Requests](https://github.com/RWTH-EBC/TEASER/pulls).
The workflow for changes is described in our [Wiki](https://github.com/RWTH-EBC/TEASER/wiki).

## How to cite TEASER

+ TEASER: an open tool for urban energy modelling of building stocks. Remmen P., Lauster M., Mans M., Fuchs M., Osterhage T., Müller D.. Journal of Building Performance Simulation, February 2017,
[pdf](http://dx.doi.org/10.1080/19401493.2017.1283539),
[bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite_jbps.bib)

### TEASER related publications

+ CityGML Import and Export for Dynamic Building Performance Simulation in Modelica. Remmen P.,
Lauster M., Mans M., Osterhage T., Müller D.. BSO16, p.329-336, September 2016,
[pdf](http://www.ibpsa.org/proceedings/BSO2016/p1047.pdf),
[bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite.bib)

+ Scalable Design-Driven Parameterization of Reduced Order Models Using Archetype Buildings with TEASER.
Lauster M., Mans M., Remmen P., Fuchs M., Müller D.. BauSIM2016, p.535-542, September 2016,
[pdf](https://www.researchgate.net/profile/Moritz_Lauster/publication/310465372_Scalable_Design-Driven_Parameterization_of_Reduced_Order_Models_using_Archetype_Buildings_with_TEASER/links/582ee96908ae004f74be1fb0.pdf?origin=publication_detail&ev=pub_int_prw_xdl&msrp=eEyK6WYemhC8wK7xkMEPRDO4obE4uxBN4-0BdBy1Ldwhy9FhCe1pXfNObJYubvC_aZN0IWDPf9uayBo3u79bsZvg3hzUoLoYRatES2ARH8c.B2cYwSICt0IOa7lD-4oAiEa_3TtrO-7k-1W9chuNQwr_VNMCpZ5ubSb-eY2D77rGUP4S6wS8m6vudUUbMlXbQQ.Cledgd1Q9fPp11nYGpcpKNhSS6bVTqAEXeMZPkiV3HsJxcVWTFj4Hr_jmLZ0MOzDxbDEZObcGiKfmTL_9k_59A)

+ Refinement of Dynamic Non-Residential Building Archetypes Using Measurement Data and Bayesian Calibration
Remmen P., Schäfer J., Müller D.. Building Simulation 2019, September 2019,
[pdf](https://www.researchgate.net/publication/337925776_Refinement_of_Dynamic_Non-Residential_Building_Archetypes_Using_Measurement_Data_and_Bayesian_Calibration)

+ Selecting statistical indices for calibrating building energy models. Vogt, M., Remmen P., Lauster M., Fuchs M. , Müller D.. Building and Environment 144, pages 94-107, October 2018. [bibtex](https://github.com/RWTH-EBC/TEASER/tree/master/doc/cite_be.bib)

+ The [Institute of Energy Efficiency and Sustainable Building](https://www.e3d.rwth-aachen.de/go/id/iyld/?) published a parametric study of TEASER where all functions and  parameters used in TEASER are gathered and explained. The publication can be found [here](https://publications.rwth-aachen.de/record/749801/files/749801.pdf).


## License

TEASER is released by RWTH Aachen University, E.ON Energy
Research Center, Institute for Energy Efficient Buildings and Indoor Climate,
under the
[MIT License](https://github.com/RWTH-EBC/TEASER/blob/master/License.md).

## Acknowledgements

This  work  was  supported  by  the  Helmholtz  Association  under  the  Joint  Initiative  “Energy System 2050 – A Contribution of the Research Field Energy”.

Parts of TEASER have been developed within public funded projects
and with financial support by BMWi (German Federal Ministry for Economic
Affairs and Energy).

<img src="http://www.innovation-beratung-foerderung.de/INNO/Redaktion/DE/Bilder/Titelbilder/titel_foerderlogo_bmwi.jpg;jsessionid=4BD60B6CD6337CDB6DE21DC1F3D6FEC5?__blob=poster&v=2)" width="200">
-->
