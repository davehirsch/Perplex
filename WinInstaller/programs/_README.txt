This version of Perplex is for testing.  There is a good deal of new code that I wrote, and I'd like help to see if it works correctly.  The impetus for writing all this code is to alter the way Perplex creates graphics files.  As you know, Perplex writes PostScript files.  You probably also know that Adobe Illustrator is essentially the only way to view these with their text intact.  In consulation with Jamie, I have altered the graphics code to use the Cairo graphics rendering library, instead of writing Postscript directly.  Thus (for now at least), Perplex requires Cairo to be installed on your machine and for the Cairo library to be accessible to the dynamic loader on your system.  The result is that not only can you have Perplex create legal PostScript, you can now have it make PDFs or SVGs directly!  This is set by a line in the plot options file.

I need a few guinea pigs to test out the new graphics code to hunt down non-obvious bugs.

These are a set of Perplex programs I compiled and linked using gfortran Windows XP.  They should run on all Windows versions back to XP, and possibly older than that.

Note that both options files have been changed somewhat.  I've supplied new options files that will work with these programs.  Your results if you run these with previous option files are undefined (meaning: who knows what might happen?)

At this time, there are only a couple of data files included (for testing).  You are encouraged to download the most current version of the data files (these change frequently)!

Cheers,
Dave Hirsch
Western WA Univ.
