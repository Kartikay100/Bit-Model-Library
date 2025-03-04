# Overview:
This repository contains the Python implementation of a series of drill bit models taken from published literature.  A full description of the models and implementation can be found in [SPE-223715-MS](https://doi.org/10.2118/223715-MS), presented at the 2025 SPE/IADC Drilling Conference in Stavanger, Norway.

The following four models are current implemented:
1. Ziaja's model ("Mathematical Model of the Polycrystalline Diamond Bit Drilling Process and Its Practical Application", 1985)
2. Detournay and Defourny model ("A phenomenological model for the drilling action of drag bits", 1992)
3. Gerbaud et al. model ("PDC Bits: All Comes from the Cutter/Rock Interaction", 2006)
4. Che et al. mode ("Chipping and crushing mechanisms in orthogonal rock cutting", 2016)

All four models may be used to estimate rate of penetration (ROP).  Weight on bit (WOB) & torque on bit (TOB) may also be calcualted using Detournay et al. and Ziaja's models.


## About this code:
- It can be employed for predicting the ROP based on 4 different models.
- It can be utilized easily for studying and comparison purposes.


## Documentation:
- The supplied model is written in Python. An example usage included in the `Example` folder:
  - To run the code, run the `SameBitModel.ipynb` notebook.
  - `Input_File.xlsx` provides an example formatting for the inputs to be compatible with the code.
  - `Output_File.xlsx` provides an example of the models outcome.
- `Parameters.txt` provides the necessary parameters used in all models with related units to facilitate understanding the code.


### Remember to:
- The path for the input and output file should be specified in the `BitModel.py`. 


### Assumption made by the authors:
- Gerbaud et al. & Che et al. models are both single cutter models and integrated into full bit model using the number of cutters (nc)
