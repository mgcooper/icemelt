Matt,

This is the version I originally gave to Matt Hoffman in
2005.

This version has been set up to run hourly time steps.

MicroMet (running within SnowModel) can be used to provide
the distributed met forcing (after this icemelt.f code has
been coded to run in distributed mode).

It looks like this code version does not write out GrADS
.gdat files.  The .ctl files are here (the .ctl variable
descriptions are not correct), but not the data outputs
and the code does not include the data-writes.  If you
want to write out .gdat files, then I can give that code
to you.  (Glen, how I did this is here:
/data1/working/people/mernild/2008_sep_visit/icemelt/).

Glen

