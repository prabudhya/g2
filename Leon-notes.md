# g2
Leon Otis' notes

8/3/17

The code necessary for running a g2 measurement is contained in TimeIntervalCounter.m and TimeIntervalCounterFunctionPool.m. TimeIntervalCounter.m corresponds to the GUI created in TimeIntervalCounter.fig while TimeIntervalCounterFunctionPool.m contains the functions that execute as part of the callbacks when buttons are clicked on the GUI.

The code was written for and tested with a Model SR620 Time Interval Counter from Stanford Research Systems connected to the computer through an Agilent GPIB interface. When the GUI is initialized, a GPIB object is created in MATLAB to handle communication with the instrument and this code will fail if some other hardware is being used for the interface. The object can then be used to send commands to the instrument and a detailed list of allowed commands is found in the SR620 manual. 

> In practice, I’ve found it is helpful to first use MATLAB’s Instrument Control tool to search for the SR620 and connect to it to check communication is working before running the code.

The g2 measurement requires the outputs of two SPCMs to be connected to the A and B channels of the SR620. As a result, any tracking of NVs and scanning must be done with the count mode of the SR620. The code for these functions is included in modified versions of ImageNVC.m and ImageFunctionPool.m. A radio button in the ImageNVC GUI must be selected and the TimeIntervalCounter GUI must also be running in order for the tracking and scanning code to work. The original functions for scanning and tracking have been modified to check if the SR620 is being used and send appropriate commands to collect count rate data in that case. As of 8/3/17, both tracking and scanning require high laser power (~.5mW at center of 4f telescope) to obtain enough counts after losses from the beamsplitter between the SPCMs. Scanning is also considerably slower than with the DAQ.

A measurement can be run from the GUI by simply entering a number of time intervals to be collected and clicking the “Measure g2” button. The program enters a loop where it sends commands to the SR620 to put it in a binary dump mode that outputs a specified number of individual time intervals measured between channels A and B. The binary data for each batch of data points is read out, converted to decimal time intervals and stored in a text file. The loop continues until the total number of data points specified by the user has been collected. There is also a stop button for ending the loop earlier.

> Time intervals are currently collected in batches of 20,000 data points and the program tracks the NV after every 100,000 points collected. These numbers are mostly arbitrary and were chosen to avoid timeouts during initial testing with low count rates from NVs. 

> When there are not many photons, the SR620 may take too much time to collect the requested number of time intervals and commands to it will time out, leaving the instrument stuck in binary dump mode. In this situation, I have found that the SR620 will remain unresponsive to remote commands unless restarted. If there is no risk of timing out, up to about 65,000 data points can be collected with each binary dump mode command though I’ve found in testing that the number of points collected in each iteration doesn’t substantially affect the overall runtime.

Once data is collected, the GUI allows the user to specify a range of time intervals and a bin width for making a histogram using a selected data file. A small separate program g2Plotter.m has been included for the same task, but without needing the SR620 connected to the computer in order to run the GUI.

I have tried to include thorough comments in the code to ensure all the implementation details are clear. 

> The main pieces of future work I had in mind were running longer measurements of at least several million time intervals to check if a g2 signal appears and investigating whether the program can be changed to handle command timeouts without problems and maybe check if tracking the NV is necessary instead of always tracking.

