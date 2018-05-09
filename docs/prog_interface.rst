overview of external interfaces
===============================


customization
*************

Most system calls to programs can be modified at run-time using ``$HOME/xoptrc``.
The program's default values are set in ``default.f90``.

All system calls (besides for conical intersection optimizations) are handled by ``getgrad.f90``, if further
modification are needed.


custom interface
****************
The general external interface (GEI) option allows
the user to call arbitrary programs and feed the gradient into ``xopt``


ORCA
****
Requirement:
 * input file: orca.in
 * ``orca30`` startup script to run orca (provided)

It is necessary that the following input lines are included::

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
 * input file: g.in
 * run-g09 startup script for Gaussian (provided)

Mopac201X
*********
Requirement:
 * binary name: mopac2012
 * input directives using a SETUP file

SETUP must contain::

  1SCF GRAD XYZ aux(42,PRECISION=9,MOS=-99999,COMP)

use of PRECISE is recommended

gamess
******

PSI4
****

Amber (sander)
**************


conical intersection optimization
*********************************
Penalty function-based (no non-diabatic coupling) CI optimizer
following Levine/Martinez DOI: 10.1021/jp0761618.
We can do Gaussian, Turbomole(dscf,ricc2), Orca and numerical gradients.

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

