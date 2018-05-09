constraints & restraints
************************

The implementation in xopt has the following limitations: Constraints are only possible for Cartesian coordinates, while restraints can be defined for internal (primitive) coordinates, e.g. bond lengths and angles.
Of course, a very restrictive restraint on an internal coordinates leads to a practical constraint. Thus, using high values for the restraining potential, also internal coordinates can be 'constrained'.

The restraining potential is in the form of a simple, quadratic function with the restraining energy $R$ being defined as

.. math::

 R(x)= k_r\cdot(p-p_{ref}) ,

with :math:`k_r` as the force constant of the restraint, $p$ the type of internal (primitive) coordinate (bond, angle or torsion) and its reference value :math:`p_{ref}`. 

The actual value of :math:`k_r` needs to be chooses carefully and is structure and level of theory specific.
The reference value :math:`x_{ref}` is taken from the initial input geometry (or from the xopt.restart file if requested).


xopt checks for the file "xopt.control" at the startup.
Its content will overwrite any commandline options!



legend::
 i,j,k,l = <integer>, atomic numbers starting from 1
 F = <floating point> eg. "1.0" 


Restrained optimizations::

 $restr
 atom i F 
 bond i j F
 angle i j k F
 dihed i j k l F
 vol sphere [vdw,<float>]

An additional gradient from the restraining potential (:math:`G_r`) will be added to the total gradient: :math:`G_total`: :math:`G + G_r`. gnorm and all other values are calculated using
:math:`G_total` ! Thus is it better to rely on energy convergence for restrained optimizations.

One can add a target value for dihedral restraints with the keyword target::

  dihed i j k l F target T

, where :math:`T` states the target value in degree.
  
  

Cartesian constrains::

 $freeze
 Hopt            #just optimize hydrogen positions, freeze the rest
 atom i          #freeze atom number i
 elem <string>   #freeze all elements <string> (lower case !)
 list <list>     #freeze a list of integers, eg. 'list 1-3,8,9,33-34'
 region [sphere,box] [atom,space] [<int>,<float>,<float>,<float>]  # freeze a region centered at an atom or a point in Cartesian space.


For constrains the gradient components (and respective entries in the B-matrix) are set to zero.

