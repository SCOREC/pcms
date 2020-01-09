# wdmapp_coupling

This git repo consists of three pieces of software

1. gene_adios2
The base code was gotten from Gabrielle as the version used for the in-memory coupling.
Its coupling capabilities were confirmed in December 2019 and still exist on cori. 
Any modification for study on XGC-GENE coupling will be put into this repo and 
eventually merged to its git repo some day.

2. XGC_adios2
This is the XGC version that serves as the second piece of the coupling. 
It was also gotten from Gabriele Merlo and its coupling has been verified.


3. wdm_coupling/cpl
This is the pass-through coupler - executable piece of code that contains the adios2 implementation of the coupling to transfer data from and to each code.

4. wdmapp_coupling/drv
This fourth piece is the obsolete interface library between XGC and GENE. 
It contains the data transfer capabilities that were initially within GENE.


This layout is subject to change upon review.
For more information, please check https://github.com/SCOREC/wdmapp_coupling/wiki

