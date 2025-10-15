# Acknowledgements
*This repository includes code and Altium files created by Yiou He and Suzanne O'Meara as part of their respective theses at MIT under the supervision of Professor David Perrault. Their works' original files, architecture, procedures, and data are the foundation to the present implementation in "MatlabWorkingDirectory". The original Altium files & MCU code from Dr. He, and the MATLAB code transcribed to the best of my ability from the thesis PDF, are included in the OriginalWork directory.*

Below is the content directly from the thesis work:
 - The MCU_and_Altium_files directory contains both the Altium files and MCU code for the 2017 1st generation HVDC power supply and the 3rd generation SCMLI HV AC power supply.
 - The MATLAB_code directory contains (currently) the LLC and RTC Boost driver and function scripts, as well as the SCMLI weight script and approximate full excel data files for CoreSizeData.xlsx, MOSFETs and Capacitor Masses.xlsx, and CoreLossData.xlsx. **The voltage multiplier weight scripts and data are not included.**

## Details on the specific procedures of each of the scripts is given in the theses available below:
- [On the Use of Switched-Capacitor Multi-Level Inverters for Electro-Aerodynamic Thrust Applications - Suzanne O'Meara](https://dspace.mit.edu/handle/1721.1/130195)
- [Towards Lightweight High-Voltage Power Conversion - Yiou He](https://dspace.mit.edu/handle/1721.1/128589
)

# HV Power Supply Weight Optimization (MatlabWorkingDirectory)
This is a script modified from the 2020 thesis by Yiou He: "Towards High Voltage Power Conversion", that has been modified for expansion of datasets from the thesis info, as well as adding volume to the result table, fixing some issues with winding patterns and consistency across driver scripts, etc. The content was parsed with OCR from the original thesis PDF and then manually checked to remove errors, then manually edited to add or change specific sections. Comments were added to explain functions of much of the code for my own understanding.

## Current issues I have identified: (delete as needed)
- I'm not entirely sure why the relative permeability of each material is fixed, and where the numbers come from.
-  Make the RTC script match the inductor and transformer design script so that every iteration of the loop writes to a row of the final design instead of to a global variable.
-  General fixes need to be performed on the RTC script to get it working.
-  Need to ensure realism of inductor data through cross-checking with more of the thesis results. Transformer comparison has already been done and is in the comparison .xlsx file.
-  A big issue is that between results, I need to change ranges like the max layers, increment Q factor, etc. not because they are important parameters for me, but because I can't have them too high due to memory limitations and I can't have them too low to risk missing designs.

## Recent changes:
- Allowed internal automatic computing for whether litz wire or solid magnet wire is better, and chooses corresponding results for the optimal choice
- Both winding pattern options work now, center and double-leg, based on the thesis text.
- U and UR core implementations for the winding pattern and insulation were included, since they weren't before.
- Multiple core densities can be used in candidate design sweep
- Interlayer tape layout option was added.
- Datasheet accurate coresizedata was added.

## RTC Script:
- Works

## Future Work:
- Import Voltage Multiplier weight scripts
- Create a way to scrape coresize datasheets online for core geometry dimensions, or just to scrape the PDFs for the correct values.
- Input core permeability data plots vs. frequency and temperature or through a plot digitizer instead of relying on the initial permeability. (Or use inductance factor of the core, or use complex permeability vs. frequency)
- Input B-CL and B-H curve data at each frequency through plot digitization instead of interpolation with 2 points.
- Compute inductor and transformer weights as a whole in the same script at once, and optimize for their coupled weight.
- Add more core geometries and core materials.
- Add U and UR cores geometries.
- Add interlayer tape to inductor.
- I want to add a plotting option for all of the candidate designs to a few graphs instead of a table, for visualization.

## Summary

### Sweeps input parameters, core size, core material, wire build, and winding pattern.
- Electrical, thermal, and mechanical parameters of each design are calculated. 
    - Electrical: Max flux density, inductance of the inductor, magnetizing impedance of the transformer.
    - Thermal: Core loss (Steinmetz), copper loss (Dowell), absolute temperature
    - Physical: Packing factor (insulation), winding width and height, interlayer tape, winding arrangement

### Rules out designs with constraints: (With Examples)
   - Max flux density >= 0.75*BSAT
   - Max temp >= 90 degC
   - Total loss 2% for inductors, 5% for transformers
   - Overall Packing factor >= 0.7,
   - Windings must fit in window height and width
   - Weight must be less than specified max
   - Layers and windings per layer must be less than max specified
   - Airgap must be greater than minimum

### Once satisfied, weight must be calculated:
   - Core
   - Wire Copper
   - Wire insulation
   - Core insulation
   - Interlayer Tape

### Weight rules of thumb for transformer & inductor:
   - The core should account for 50% of the weight,
   with the other 50% taken up by the copper.
   - Current density of the secondary is usually smaller
   than the rule-of-thumb 500A/cm^3.
   - Core loss %: core loss % corresponds to less current
   density in the wire

## Common Errors

1. Ensure the excel files are not open elsewhere, or else you will get an error: 
"Unable to open file 'path' as a workbook. Check that the file exists,..."

2. If your variables are too restrictive for your data, you will get a table with
no values, as no set of variables match your requirements

3. When running the script, if the weight seems implausibly high, check the results table for both inductor and transformer, and observe if any of the values in the table approach the bounds of the parameter ranges. If they do, expand those ranges. Often, expanding those ranges allows for more optimized cores. Once each of the parameters is in a good band within its range, the optimal design can be certain to have been achieved.

## Data formatting:

You need to collect geometry and size data, as well as material-specific data.

### SizeData:

Parameter meanings:

 PriW = primary winding width
 PriH = primary winding height
 SecW = secondary winding width
 SecH = secondary winding height

Dimensional references:
 Lc = center-leg width
 T = core thickness (used for EE or U cores)
 rAc = equivalent core radius (used for ER or UR cores)

For geometry: You need geometries with known dimensions
You need:
- Core effective length Le (mm)
- Core window area height and width (mm^2)
- Core volume (mm^3) as Le*Ac
- Center leg and Side leg area

All units are mm, or mm^2, or mm^3

- Columns go as follows:

rowNumber, coreName, Ve, Ae, Le, coreShapeIndex, empty, priWindingWidth, priWindingHeight, secWindingWidth, secWindingHeight, coreWindowWidth, coreWindowHeight

### LossData:

For material: You need flux density and power loss at discrete known frequencies, saturation flux density, and permittivity.
Freq is in Hz, Bfield is in T, Ploss is in mW/cm^3, BSAT is in T, MU is in H/m

- Freq: at least one frequency must be within 40% of operating freq. Frequency is in pairs of values. cols in pairs of 2 freq. per material
- Bfield: exactly 2 B points per freq. This is by design, so choose 2 points close to your intended operating point or just in the middle of the CL loss data to avoid the low and high ranges, and it should be more accurate than any regression fit with more points (in theory)
- Ploss: 2 loss values corresponding to B values, in mW/cm3
- BSAT: design derates to 0.75 of value. Put one value in C per material
- MU: material mu_r in column C

Core loss data in CL provided in manufacturer datasheets 



