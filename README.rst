CIAO: Calern Imaging Adaptive Observatory
=========================================

A couple different pieces of software used to run the prototype AO
system called CIAO, being developed for the C2PU-Epsilon telescope.

Organisation:
-------------

The software architecture is built around the idea of shared memory
data structures, identical to the ones created in the context of the
SCExAO project at the Subaru Telescope.

Any AO system revolves around three components:

- A deformable mirror to control the wavefront the DM)
- A camera used to sense the wavefront (the WFS)
- A supervisor that looks at the WFS and drives the DM

Consequently, the software architecture follows this split:

- A server runs to control the acquisition by the camera
- A server runs to control the commands sent to the DM
- A main GUI and a few auxilliary tools are used to control everything

And an additional directory containing libraries that are currently
not installed system-wide, either because too specific or not
incredibly important.

