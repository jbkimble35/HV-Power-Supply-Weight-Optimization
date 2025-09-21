# HV-Power-Supply-Weight-Optimization
This is a script modified from the 2020 thesis by Yiou He: "Towards High Voltage Power Conversion", that has been modified for readability and easier data importing, as well as taking into account Volume, fixing some issues with winding patterns and consistency across driver scripts, and more to come. The content was parsed with OCR from the original thesis PDF and then manually checked to remove errors. Comments were added to explain functions of much of the code that I deemed confusing.

Current issues I have identified: (delete as needed)
-  Check window area calculations and if it's being doubled for the 2 window cores and vice versa. I believe this is an issue case because of the data mismatch in my excel file.
-  Make the RTC script match the inductor and transformer design script so that every iteration of the loop writes to a row of the final design instead of to a global variable.
-  The CoreLossData.xlsx and CoreSizeData.xslx are too large to run on my home computer, the n-dimensional array generated takes too much RAM.
-  Need to change the winding arrangement code according to thesis 147-158 to allow for ER double-leg as well as U and UR geometries (since EE,ER,U,UR are the best for HV power magnetics)
-  Need to implement winding pattern == 2 code for the transformer script.
-  General fixes need to be performed on the RTC script to get it working.
-  Main LLC script doesn't work for high voltages for the inductor, which is the main purpose of the script. I believe this is due to some airgap or geometry issue.
-  Need to ensure realism of data through cross-checking with more of the thesis results. Transformer comparison has already been done and is in the comparison .xlsx file
-  Need to comment what each block of code is doing for readability and debugging
-  Need to add error handling and debugging printing in the terminal for common issues and diagnosis of problems.

Currently working on:
- I am going to edit this script so that it takes into account higher-voltage transformer design parameters, such as volts-per-layer, epoxy potting, and interlayer tape.
