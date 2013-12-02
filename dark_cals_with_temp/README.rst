This contains a hacked version of the dark cal warm pixel correction
script (which is warm_pix_corr.py in the aca_dark_cal project).  The
hacked version puts out a temperature for all dark cals and writes
everything out to a text file 'dark_cals_with_temp.txt'.

The temperatures are from the mica HDR3 data when available (or when
it worked with the starcheck-reading script to determine replica
times) and AACCCDPT from the engineering archive when HDR3 data could
not be retrieved.  The field 'temp_source' in the text file notes
which temperature data was used.  

I've tagged the 1999:223 data with "-10" and 'GUESS" as I did not have
eng archive temps for that either.
