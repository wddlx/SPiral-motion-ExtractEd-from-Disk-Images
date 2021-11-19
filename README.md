# Spiral-arm-motion-calculator
This code aims to offer an easy way to calculate the location of spiral arms of a protoplanetary disk and calculate the motion between different epochs.

The main code is the file 'calculator.py', please import it while using. 
The 'example.py' and the 'example.ipynb' are examples to deproject a disk figure, fit the locations of the centers of spiral arms, and calculate the motion between figures from different epochs.

Your input width and regions will decide what you get, so try to use "Ispt" function to check the fitting.

If you run into an error "'ydata' must not be empty", check if the center+width larger than y_fin or center-width smaller than y_init. Since the default center is the brightest location, make sure that in the selected region, the center of the spiral arm is also the brightest. e.g. Don't include brighter regions like the inner ring. You may use "Ispt" function to check it. 

I'm still developing the package, contact me if you run into any problem. 


Thanks a lot to @Bin Ren and @bobpanda for the help.