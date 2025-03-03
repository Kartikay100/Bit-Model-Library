*Overview:* 
This Code main function is to calculte the ROP from the 4 following models:
1- Ziaja's model ("Mathematical Model of the Polycrystalline Diamond Bit Drilling Process and Its Practical Application", 1985)
2- Detournay and Defourny model ("A phenomenological model for the drilling action of drag bits", 1992)
3- Gerbaud et al. model ("PDC Bits: All Comes from the Cutter/Rock Interaction", 2006)
4- Che et al. mode ("Chipping and crushing mechanisms in orthogonal rock cutting", 2016)

It is also capable of calculating the WOB & TOB using Detournay et al. and Ziaja's models.


*About this code:*
- It can be employed for predicting the ROP based on 4 different models.
- It can be utilized easily for studying and comparison purposes.


*Documentation:*
- The supplied model is written in Python. To run the code, run XXXXX.
- Attached with the code two files which are Input_File.xlsx and Output_File.xlsx.
   . The former is showing how the inputs should be fed in to be compatible with the code.
   . The later is a sample of the models outcome.
- Attached also Parameters.txt, illustrating all the parameters used in all models with related units to facilitate understanding the code.
- Paper pre-print for test cases and validation.


*Remember to:*
- Do not forget to modify the path of this files at the end of the code.


*Assumption made by the authors:*
Gerbaud et al. & Che et al. models are both single cutter models and integrated into full bit model using the number of cutters (nc)