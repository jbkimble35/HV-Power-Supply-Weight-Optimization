# Acknowledgements
*This repository includes code and Altium files created by Yiou He and Suzanne O'Meara as part of their thesis work at MIT under the supervision of Professor David Perrault. Their work's original files, architecture, procedures, and data are the foundation to the current implementation. To the best of my ability in transcribing the code from PDFs, the unedited code and original Altium files are included in the OriginalWork directory.*

## Details on the specific procedures of each of the scripts is given in the theses available below:
- [On the Use of Switched-Capacitor Multi-Level Inverters for Electro-Aerodynamic Thrust Applications](https://dspace.mit.edu/handle/1721.1/130195)
- [Towards Lightweight High-Voltage Power Conversion](https://dspace.mit.edu/handle/1721.1/128589
)

# HV Power Supply Weight Optimization
This is a script modified from the 2020 thesis by Yiou He: "Towards High Voltage Power Conversion", that has been modified for readability and easier data importing, as well as taking into account Volume, fixing some issues with winding patterns and consistency across driver scripts, etc. The content was parsed with OCR from the original thesis PDF and then manually checked to remove errors. Comments were added to explain functions of much of the code for my own understanding.

## Current issues I have identified: (delete as needed)
-  Make the RTC script match the inductor and transformer design script so that every iteration of the loop writes to a row of the final design instead of to a global variable.
-  General fixes need to be performed on the RTC script to get it working.
-  Main LLC script doesn't work for high voltages and powers for the inductor. I believe this is due to some airgap or geometry issue.
-  Need to ensure realism of inductor data through cross-checking with more of the thesis results. Transformer comparison has already been done and is in the comparison .xlsx file.
-  Need to comment what each block of code is doing for readability and debugging
-  Need to add error handling and debugging printing in the terminal for common issues and diagnosis of problems.
-  Need to check Ns4 to make sure it's not vital for operation of the script, since it's not being used in the function that it calls on the GPU.
-  Need to replace core insulation thickness with a new value, need to take into account the epoxy between wires in each layer and between layers, need to include secondary layer stress and interlayer tape, need to filter designs by their volts-per-layer.
- I need to edit the inductor script to solve some of the same issues I solved in the transformer one, i.e. maybe adding a GPU option and solid core vs. litz wire option

## Recent changes:
- Added GPU computing options for this transformer script, but the old CPU options also work
- Allowed internal automatic computing for whether litz wire or solid magnet wire is better, and chooses corresponding results for the optimal choice
- Both winding pattern options work now, center and double-leg, based on the thesis text.
- U and UR core implementations for the winding pattern and insulation were included, since they weren't before.
- Made it more efficient by allowing GPU computation on some of the most computation-intensive array functions, which were tested with the code profiler toolbox
- Multiple core densities can be used in candidate design sweep

## RTC Script:
- Currently in a transition state, doesn't work
- Currently am fixing up the RTC script so that the excel writing and global variables are removed from the parfor loop, so the RTC script largely matches the inverter transformer one. I already made the Design(1 to 33).data into Design( : , 1 to 32), erased global variables, erased some of the excel writing and runNumber stuff, and need to match the two scripts now one to one with the loop parts. I also made the driver script have a parfor loop to match the other one.

## Future Work:
- Import Voltage Multiplier weight script


