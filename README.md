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

Recent changes:
- Added GPU computing options for this transformer script, but the old CPU options also work
- Allowed internal automatic computing for whether litz wire or solid magnet wire is better, and chooses corresponding results for the optimal choice
- Both winding pattern options work now, center and double-leg, based on the thesis text.
- U and UR core implementations for the winding pattern and insulation were included, since they weren't before.
- Made it more efficient by allowing GPU computation on some of the most computation-intensive array functions, which were tested with the code profiler toolbox
