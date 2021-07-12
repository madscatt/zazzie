Add the following lines to your .vmdrc to make it easy to load the packages; replace "/home/paul/scripts/" with the directory in which you unpacked the packages:
  lappend auto_path /home/paul/scripts/la1.0
  lappend auto_path /home/paul/scripts/orient
Aligning a molecule to its principal axes

Load the molecule into vmd, and run the following commands to align the first, second, and third principal axes to the x, y, and z axes.
  package require Orient
  namespace import Orient::orient

  set sel [atomselect top "all"]
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 2] {0 0 1}]
  $sel move $A
  set I [draw principalaxes $sel]
  set A [orient $sel [lindex $I 1] {0 1 0}]
  $sel move $A
  set I [draw principalaxes $sel]
You may want to use the principal axes (stored in I) for other purposes than alignment, or modify the script to get information about the moments of inertia.


