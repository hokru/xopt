parallel
========

Parallel numerical gradients (max. 3*atoms threads) are available through an
add-on programm called ``xopt.pgrad`` (see subdirectory \textit{pgrad}).
It can be compiled by running ``make pgrad` or by running ``make`` in the
pgrad directory itself. It has only been tested with OpenMPI.

Communication with xopt is handled through 2 files:
``xopt.para.tmp`` contains the input for ``xopt.pgrad``.
It hold the number of atoms, the xyz coordinates (\AA\ and the atom type (as integer).
For water it reads::

 3
 -1.769142     -0.076181      0.000000 8
 -2.065745      0.837492      0.000000 1
 -0.809034      0.001317      0.000000 1


``xopt.pgrad.tmp`` contains the output from ``xopt.pgrad``.
Two gradient matrices are written after each other.
There are two gradients since in the case of conical section (CIopt) optimizations one
can read the energy from two states (eg. GS and S1) from the same output and form both
numerical gradients simultanously.
If there are not two energy values inside xopt.energy.tmp, then the 2nd gradient
in ``xopt.pgrad.tmp`` will be zero. 
The gradient is in atomic units and the file is organized as follows::

 grad1(x,atom 1), grad1(y,atom 1),grad1(z,atom 1)
 .
 .
 grad1(x,atom n), grad1(y,atom n),grad1(z,atom n)
 grad2(x,atom 1), grad2(y,atom 1),grad2(z,atom 1)
 .
 .
 grad2(x,atom n), grad2(y,atom n),grad2(z,atom n)

Each thread write its output into xopt.slave.* files, enumerated with process rank.
The master process hold the rank 0. Its output will be copied to xopt.pgrad.out afterwards.


Internally \xopt will call \texttt{xopt.pgrad} as::

 mpiexec -by-node -n <nproc> -output-filename xopt.slave $HOME/bin/xopt.pgrad 

meaning the ``xopt.pgrad`` binary has to reside in the users ``$HOME/bin``. The mpiexec call can 
be customized in ``.xoptrc``.

task distribution
*****************
The parallel task assigment is kept very simple: each thread gets
``3*atoms/nproc`` tasks.

Each tasks consists of a single gradient components (e.g. 2 displacements for x-coordinate of atom 1).
This likely leaves several leftover tasks (check modulus).
These will be assigned to the process with the highest rank! 
For best parallel performance you should aim to minimize the number of leftover tasks.

In the case of 19 atoms, we have 57 tasks to calculate.
If we use 16 MPI threads, then 15 threads will calculate 3 tasks,
while the 16th thread will calculate 12 tasks (3+9 leftovers),
resulting into a severe bottleneck.
In this case it would be better to use only 14 MPI threads,
since then 13 threads will run with 4 tasks and the 14th with 5 (4+1 leftover).

tip: quickly check with the python interpreter for leftover tasks: (57\%16=9 and 57\%14=1).

scratch directories for xopt.pgrad
***********************************
There are possible modes::

 1.[default] xopt.pgrad will automatically make scratch directories in the working
 directory named \verb|xopt-tmp-$(PID)|

 2.[set by -scratch ] You can set a custom scratch directory (like /scratch/myname/) for computations on clusters
 by adding an additional line to \verb|xopt.pgrad.tmp|. The above water example would then read:

 3
 -1.769142     -0.076181      0.000000 8
 -2.065745      0.837492      0.000000 1
 -0.809034      0.001317      0.000000 1
 /scratch/username


Xopt can take care of all that by itself.
The first mode is the default and used when only ``-numgrad -n <nproc>`` is specified.
The second mode is actived if additionally ``-scratch ...`` is set.
Additionally, xopt will do the energy computation in the same scratch path.
This way, you can submit the computation directly on the NFS directory on the cluster and the actual
computations will be carried out in the specified local scratch directories.

As long as ``mygrad.sh`` and ``xopt.pgrad.tmp`` are present, you can call ``xopt.pgrad`` directly (with mpirun...)
if you want to use it in your own scripts without xopt.

