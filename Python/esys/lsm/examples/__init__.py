#############################################################
##                                                         ##
## Copyright (c) 2003-2017 by The University of Queensland ##
## Centre for Geoscience Computing                         ##
## http://earth.uq.edu.au/centre-geoscience-computing      ##
##                                                         ##
## Primary Business: Brisbane, Queensland, Australia       ##
## Licensed under the Open Software License version 3.0    ##
## http://www.apache.org/licenses/LICENSE-2.0              ##
##                                                         ##
#############################################################
"""
Example LSM python modules and packages.

Example Simulations
===================
    MPI and C{mpirun}
    -----------------
        The Lattice Solid Model implementation takes advantage of
        the MPI (Message Passing Interface U{http://www.mpi-forum.org/})
        to provide parallel execution. MPI programs
        are executed in an MPI environment, usually by invoking the C{mpirun}
        command.

        Unfortunately, different implementations of MPI, (eg LAM, SGI's MPT,
        MPICH, etc) do not necessarily share the same command-line options
        for the C{mpirun} command. To make running LSM simulations a
        little easier for users, the ESyS-Particle software provides
        some wrappers scripts which invoke the appropriate C{mpirun}
        command.

    Running Simulations on the Altix 3700
    -------------------------------------
        The folling sub-sections explain the details of running
        LSM simulations on the ess.esscc.uq.edu.au,
        an SGI Altix-3700 supercluster. In the examples, running
        the C{esys.lsm.examples.waveprop.WaveSim}
        particle-wave-propagation python-module is used to illustrate
        command syntax.

        Initialising the user environment - environment modules
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            On the Altix 3700, ess.esscc.uq.edu.au, I{Environment Modules}
            U{http://modules.sourceforge.net/} is used for setting up
            user environments for various locally installed software.
            To define the assorted C{module} commands and make the various
            packages available, at a bash prompt execute::
              user@ess$ source /opt/modules/default/init/bash
              user@ess$ export MODULEPATH=$MODULEPATH:/raid2/tools/modulefiles
    
            Generally, it is a good idea to put these two commands into a
            C{~/.bashrc} file as follows::
              #
              # Initialise the environment-modules stuff
              #
              if [ -f /opt/modules/default/init/bash ]; then
                  source /opt/modules/default/init/bash
                  export MODULEPATH=$MODULEPATH:/raid2/tools/modulefiles
    
                  #
                  # Load environments for various favourite modules.
                  #
                  module load gnuplot-4.0.0
                  module load pbspro5.4
              fi
    
            To set up the environment for running LSM simulations, execute::
              user@ess$ module load lsm/sgimpi
    
            If the environment has been successfully loaded, executing::
              user@ess$ which mptrunLsm
    
            should return a valid path to the C{mptrunLsm} script.
        
        Running esys.lsm.examples.waveprop.WaveSim - mptrunLsm
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            To run the L{esys.lsm.examples.waveprop.WaveSim} python module
            using the SGI's MPT (Message Passing Toolkit), execute the
            following at the command prompt::
              user@ess$ mkdir MyWaveSim
              user@ess$ cd MyWaveSim
              user@ess$ nice -n 15 mptrunLsm esys.lsm.examples.waveprop.WaveSim
    
            This command will begin executing a wave propagation simulation,
            generating some output to stdout and creating various data files
            in the current working directory.
            Notice the use of C{nice}, to lower the execution priority
            of the simulation, this is a courtesy to other users when running
            longer multi-cpu jobs in the interactive CPU set on multi-user
            systems.
            The C{mptrunLsm} script simply examines the first command line
            argument (treating it as a python module name) and B{exec}'s the
            associated python-script
            ($LSMDIR/lib/python2.3/site-packages/esys/lsm/examples/waveprop/WaveSim.py)
            as an MPI python process with appropriate C{mpirun} command line
            arguments.
    
            The L{esys.lsm.examples.waveprop.WaveSim} module accepts several
            command line options which control various parameters of the
            simulation. The C{--help} option::
              user@ess$ mptrunLsm esys.lsm.examples.waveprop.WaveSim --help
    
            causes the listing of the available command line options including
            descriptions of options and default values.

        Visualisation
        ~~~~~~~~~~~~~
            The wave-propagation simulation writes particle position and
            displacement data to C{particle_*.txt} text files. Two-dimensional
            wave propagation I{snap-shots} can be visualised using the
            L{esys.lsm.examples.DisplacementPlotter} module. To visualise
            the displacement data in file C{particle_9.txt}, execute the
            following at the command prompt::
              user@ess$ lsmpython esys.lsm.examples.DisplacementPlotter particle_9.txt
            
            This will display a VTK (Visualisation ToolKit) window with the
            particle displacement magnitude plotted as a 3D surface. The window
            is interactive: spin the camera by dragging with left mouse button,
            translate camera by dragging with middle mouse-button and zoom the
            camera by dragging the right mouse-button and close the window by
            typing 'q'.
    
            Here, the C{lsmpython} script is used to run the
            C{DisplacementPlotter} python module.
            This module does not rely on an MPI implementation so can be
            run directly under python. The C{lsmpython} script is provided as a
            convenience for running modules specified by module-name
            (C{esys.lsm.examples.DisplacementPlotter}) as opposed to being
            specified by module-script-filename
            ($LSMDIR/lib/python2.3/site-packages/esys/lsm/examplesDisplacementPlotter.py)
            E{-} less typing and running a python-module is independent of the
            install location.
            
            The L{esys.lsm.examples.DisplacementPlotter} accepts command line
            options::
              user@ess$ lsmpython esys.lsm.examples.DisplacementPlotter --help
    
            will print valid command line options, descriptions and default
            values.
    
            A L{esys.lsm.examples.waveprop.WaveSim} simulation also produces
            seismograph data files. This seismograph data can be visualised
            with C{gnuplot} U{http://www.gnuplot.info/}. For example, running
            the WaveSim module with default parameters produces the seismograph
            data file C{srcToCorner_265.000_229.631_0.000.txt}.
            This file contains data for a seismograph placed at coordinate
            M{(265.000,229.631,0.000)} (this seismograph is a member of a
            group of seismographs which are positioned on a line extending
            from the source disturbance to a corner of the particle block,
            hence the "srcToCorner_" file name prefix).
            Each line of the file contains 10 values separated by whitespace:
            time, displacement (dx,dy,dz), velocity (vx,vy,vz) and
            acceleration (ax,ay,az).
            The M{x} component of the displacement can be plotted using::
              user@ess$ module load gnuplot-4.0.0
              user@ess$ gnuplot
                      .
                      .
                      .
              Terminal type set to 'x11'
              gnuplot> plot 'srcToCorner_265.000_229.631_0.000.txt' using 1:2 with line
              
            which produces a line plot of M{x}-displacement (column 2 of data
            file) versus time (column 1 of the data file). To gnuplot the
            displacement magnitude M{sqrt(S{delta}x*S{delta}x +
            S{delta}y*S{delta}y+S{delta}z*S{delta}z)}::
              gnuplot> plot 'srcToCorner_265.000_229.631_0.000.txt' using 1:(sqrt($2*$2+$3*$3+$4*$4)) with line

            Another dataset produced during the simulation is "record-section" data.
            Again, this is seismograph displacement, velocity and acceleration
            data, but it is for multiple seismographs. The M{x}-displacement
            record-section data can be displayed in C{gnuplot} using::
              gnuplot> plot 'srcToCorner_record_section_reordered.txt' using 1:($2*0.0004 + $3) with line

            This plots multiple seismograph datasets on a single set of axes and
            allows the visualisation of P-wave and S-wave propagation.
            The 'srcToCorner_record_section_reordered.txt' file contains 11 columns
            of data (time, distance-to-source, displacement(dx,dy,dz),
            velocity(vx,vy,vz), acceleration(ax,ay,az)).

        Submitting Simulation Jobs using PBS Pro
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            To set up the environment for PBS Pro U{http://www.pbspro.com},
            load the C{pbspro5.4} module::
              module load pbspro5.4
    
            To submit a script through PBS::
              qsub qsubScript.sh
              
            An example C{qsubScript.sh} script for for submitting a wave
            propagation simulation to PBS::
              #!/bin/bash
              ### Job name must be <= 15 characters
              #PBS -N MyWaveSim
              ### Specify a PBS queue, q1 - 2 to 4 CPU jobs, q2 - 8 to 32 CPU jobs
              #PBS -q q2
              ### Specify required number of CPUs
              #PBS -l ncpus=18
              ### Redirecting Output and Error files
              #PBS -e /raid2/user/MyWaveSim/out.err
              #PBS -o /raid2/user/MyWaveSim/out.err
              ### Mailing ALerts: a -abort, b -begin execution and e -end execution
              #PBS -m abe
              #PBS -M user@uq.edu.au
              
              #
              # Setup user paths to run LSM stuff
              #
              source /opt/modules/default/init/bash
              export MODULEPATH=$MODULEPATH:/raid2/tools/modulefiles
              module load lsm/sgimpi
              export LSMPYTHONPKGDIR=$LSMDIR/lib/python2.3/site-packages/esys/lsm
              
              #
              # Where we want data to be saved, create $OUTPUTDIR
              # if it doesn't exist.
              #
              export OUTPUTDIR=/raid2/user/MyWaveSim
    
              if [ -d $OUTPUTDIR ]; then
                echo Output dir $OUTPUTDIR exists.
              else
                echo Creating $OUTPUTDIR ...
                mkdir $OUTPUTDIR
                cd $OUTPUTDIR
              fi
              cd $OUTPUTDIR
              #
              # Remove any existing data files.
              #
              rm srcTo*.txt surf*.txt particle_*.txt out.err  out.log
              
              #
              # Use mpirun and dplace to place MPI processes on CPU's.
              #
              mpirun -up 128 -np 1                                  \\
                `which dplace` -e -c0-31                            \\
                `which mpipython`                                   \\
                $LSMPYTHONPKGDIR/examples/waveprop/WaveSim.py       \\
                -b[512,320,1]                                       \\
                -m8000                                              \\
                -n16                                                \\
                -s0.35                                              \\
                -d[4,4,1]

"""
from .simple import *
from .waveprop import *
from . import Wave2d
from esys.lsm.util import InstallInfo
if (InstallInfo.haveVtk()):
    from . import DisplacementPlotter
