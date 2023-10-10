# Script_repository
  Python script repository for neuroimaging analysis
Get it with `pip install git+ssh://git@github.com/d-van-de-velden/dv_code.git`



## Introduction
  First of all, I would like to highlight the three pillars on which this toolbox is built.
- Collaboration
- Simplicity
- Reproducibility

### On Collaboration:
  The functionality, and therefore the power, of this toolbox depends on the **collaboration** of a wide range of people. If you want to collaborate, click [here](https://www.skeidelab.com/).
  
### On Simplicity:
  Despite knowledge of what tools, steps and analysis you want to perform the *only* thing you have to do is create a `main_analysis.py` script to code (serves also as documentation) your computational analysis. Functionalities of this toolbox will provide funtions to use that allow you to customize your analysis flow and takes care about everything else in the background. This aims at highest **simplicity** in usage as possible.
  
### On Reproducibility:
  Due to current discussion in neuroimaging community, this toolbox adheres to the [BIDS](https://bids.neuroimaging.io/) standard when it comes to direcotry and file structure. This ensures **reproducibility** when it comes to datastructure and function appliaction.
<sub>However, important results will be provided in a detached folder aside  the BIDS structure to show your analysis history in a clear and informative way.</sub>


### What this toolbox can do
  Furthermore, the functionalities of this toolbox should be explained:
1. Downloading DICOM files and organizing a BIDS structure.
2. Preprocessing of structural and functional MR datasets.
3. Analyzing behavioral task data (D-Prime analysis, hDD-Modelling)
4. Applying approaches of 1st and 2nd level analysis on functional MR datasets. 


## Getting started
Let us construct a `main_analysis.py`. First thing you want to run after using:
`from dv_code.scripts.initialize_study import initialize_study`
`initialize_study()`-function.
