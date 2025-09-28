# HV-Power-Supply-Weight-Optimization
This is a script modified from the 2020 thesis by Yiou He: "Towards High Voltage Power Conversion", that has been modified for readability and easier data importing, as well as taking into account Volume, fixing some issues with winding patterns and consistency across driver scripts, and more to come. The content was parsed with OCR from the original thesis PDF and then manually checked to remove errors. Comments were added to explain functions of much of the code that I deemed confusing.

Current issues I have identified: (delete as needed)
-  Make the RTC script match the inductor and transformer design script so that every iteration of the loop writes to a row of the final design instead of to a global variable.
-  General fixes need to be performed on the RTC script to get it working.
-  Main LLC script doesn't work for high voltages for the inductor, which is the main purpose of the script. I believe this is due to some airgap or geometry issue.
-  Need to ensure realism of data through cross-checking with more of the thesis results. Transformer comparison has already been done and is in the comparison .xlsx file
-  Need to comment what each block of code is doing for readability and debugging
-  Need to add error handling and debugging printing in the terminal for common issues and diagnosis of problems.
-  Need to check Ns4 to make sure it's not vital for operation of the script, since it's not being used in the function that it calls on the GPU.
-  Need to include epoxy potting considerations and interlayer tape thickness additions to the window area

Currently working on:
- I am going to edit this script so that it takes into account higher-voltage transformer design parameters, such as volts-per-layer, epoxy potting, and interlayer tape.
- I need to edit the inductor script to solve some of the same issues I solved in the transformer one, i.e. adding a GPU option and solid core vs. litz wire option

Recent changes:
- Added GPU computing options for this transformer script, but the old CPU options also work
- Allowed internal automatic computing for whether litz wire or solid magnet wire is better, and chooses corresponding results for the optimal choice
- Both winding pattern options work now, center and double-leg, based on the thesis text.
- U and UR core implementations for the winding pattern and insulation were included, since they weren't before.
- Made it more efficient by allowing GPU computation on some of the most computation-intensive array functions, which were tested with the code profiler toolbox
- Multiple core densities can be used in candidate design sweep

RTC Script:
- Currently am fixing up the RTC script so that the excel writing and global variables are removed from the parfor loop, so the RTC script largely matches the inverter transformer one. I already made the Design(1 to 33).data into Design( : , 1 to 32), erased global variables, erased some of the excel writing and runNumber stuff, and need to match the two scripts now one to one with the loop parts. I also made the driver script have a parfor loop to match the other one.
