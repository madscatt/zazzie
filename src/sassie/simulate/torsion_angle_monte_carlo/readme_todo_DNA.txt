I went through and verified these two gui mimic files.  There are two issues for you to look into
I do not think I am calling the minimization correctly (look at line 667 of double_stranded_nucleic.py).  This is where the both gui mimic files break when setting `psf_flag = True` (I have it set to False in the attached files so they would run.
Every time I have run the nucleosome, monte_carlo has determined there was a collision between the flexible DNA and the rest of the molecule.  I mentioned that I forced it to run by turning off that check in monte_carlo.py, not ideal.


