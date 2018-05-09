overview of external interfaces
===============================


customization
-------------

Most system calls to programs can be modified at run-time using ``$HOME/xoptrc``.
The program's default values are set in ``default.f90``.

All system calls (besides for conical intersection optimizations) are handled by ``getgrad.f90``, if further
modification are needed.


setting custom system calls via ``$HOME/.xoptrc``

::

   call_orca = orca_scipt_or_binary
   call_gaus = run-g09
   call_psi4 =psi4
   call_gei=mygrad.sh
   call_mpi=mpiexec -bynode
   call_mopac=mopac16
   scratch=/local_scratch/
   maxiter=500  # max iterations
   maxd=0.35    # max displacement bohr
   gconv=1e-3   # gradient norm threshold
   econv=1e-7   # energy change threshold
   maxgrad=1e-3 # max. gradient threshold
  




numerical gradients
*******************
Arbitray numerical gradients can be computed via a script/executable named  ``mygrad.sh`` located in the working directory.
It does not need to be bash script. It will be excuted by a system call inside xopt.
It needs to do 3 things::

- Have the program to be executed use ``xopt.xyz`` as input
- run the program
- provide the energy output as a single number inside ``xopt.energy.tmp``


An example for Turbomole would look like::

 #!/bin/bash
 #1. make use of xopt.xyz
 babel -ixyz xopt.xyz -otmol coord 
 #2. execute
 dscf > dscf.out
 #3. write energy into file
 grep "|  total energy" dscf.out | awk '{print $5}' > xopt.energy.tmp 


``xopt`` will read the ``xopt.energy.tmp`` value for a 2-point (central) finite difference computation
of the gradient with a step size of 0.005 bohr.

general external interface
**************************
The general external interface (GEI) option allows
the user to call arbitrary programs and feed the gradient into ``xopt``. Is uses similar
technology as the above numerical gradient.

Use any script that transforms a ``xopt.xyz`` file into an energy and gradient file with the 
following format in a.u.::

 energy
 grad1(x,atom 1), grad1(y,atom 1),grad1(z,atom 1)
 .
 .
 grad1(x,atom n), grad1(y,atom n),grad1(z,atom n)



The default name of the script is ``mygrad.sh`` and the default gradient file is ``xopt.grad``




ORCA
****
Requirement:
 * input file: ``orca.in``
 * ``orca40`` startup script to run orca (provided)

It is necessary that the following input lines are included

::

   ! ENGRAD
   *xyzfile <chrg> <mult> xopt.xyz


Turbomole
*********
Prepare a standard input using define. Supported are HF/DFT calculations using ridft and dscf, and calculations using the ricc2 module (eg. for MP2.).
A 'coord' file is automatically written and updated. If ``$rij`` is found in the control file, ridft will be called instead of dscf.

Requirement:
 * ridft/dscf/ricc2 in ``$PATH``

Gaussian09
**********
Requirement:
 * input file: ``g.in``
 * ``run-g09`` startup script for Gaussian (provided)

Mopac201X
*********
Requirement:
 * binary name: ``mopac2016``
 * input directives using a ``SETUP`` file

SETUP must contain at least:

::

  1SCF GRAD XYZ aux(42,PRECISION=9,MOS=-99999,COMP)

use of ``PRECISE`` is recommended.

gamess
******

`experimental`

nwchem
******

`experimental`

PSI4
****

needs to be adapted to changes in psi4 v1.2

Amber (sander)
**************

requires modified sander


conical intersection optimization
*********************************
Penalty function-based (no non-diabatic coupling) CI optimizer
following Levine/Martinez DOI: 10.1021/jp0761618.
We can do ``Gaussian``, ``Turbomole(dscf,ricc2)``, ``Orca`` and (modified) ``Amber`` and numerical gradients.

You need to make 2 directories named ``stateJ.xopt` and ``stateI.xopt``.
Prepare the input for each state inside the directories.
It should work for Turbomole, ORCA and G09 if you follow the general preparation guidelines above.

We assume state I < state J, e.g. J=I+1.


The ``xopt`` output will print something like:

``gap[eV]:   0.024   penalty:  13.1 E(low):   -546.9431436 E(high):   -546.9422645 root flip: F``,
where E(low) denotes the lower state (eg. groundstate) calculated as state I in stateI.xopt and E(high)
as the higher state calculated as state J inside stateI.xopt.

If at any stage during the optimization E(low)>E(high) (I>J)
root flip will be set true (=T) and E(low/higher) will be interchange, e.g.
E(low) will be the energy obtained in ``stateJ.xopt``.
This is checked for each optimization step individually and is not tracked through previous steps.

The strategy to increase the penaly :math:`\sigma` is as follows:

.. math::

   \sigma= \sigma+(2\Delta E_{ij}/\alpha) \ ,

where :math:`\Delta E_{ij}` is the energy gap between state I and J, $\alpha$ the smoothing factor (see paper).
$\sigma$ is increased when :math:`\Delta E_{ij}` is larger than 1e-3 and the penalty function change :math:`\Delta F_{ij}` smaller than 5e-5.

Note 1: for gaussian groundstate SA-CASSCF calculations add IOp(10/97=100) and IOp(5/97=100).
It switches the CI vectors internally so one gets the proper groundstate gradient
(e.g. for CASSCF(2,2,nroot=2,stateaverage)). Not sure about nroot > 2.

