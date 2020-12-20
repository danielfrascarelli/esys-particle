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
Defines classes for parsing command line options of wave propagation
simulations.
"""
import optparse
import re
import copy

def getListFromString(opt, str, mapCallable):
    """
    Parses a given string to extract a list. The expected format of the
    string is '[e1,e2,e3,...]'. Returns C{map(mapCallable, [e1,e2,e3,...])}.
    @type opt: string
    @param opt: Command line option, only used for for creating
               messages when raising OptionValueException
    @type str: string
    @param str: The option value string
    @type mapCallable: callable
    @param mapCallable: Callable object which converts a string
                       list element into another type, eg int, float, etc

    @rtype: list
    @return: list of converted-type elements.
    """
    regex = re.compile("\[(.*)\]")
    match = regex.match(str)
    if (match != None):
        try:
            return list(map(mapCallable, str.split(match.group(1), ",")))
        except ValueError:
            raise \
                optparse.OptionValueError(
                    (
                      "option {0!s}: Could not convert string elements of"+
                      " {1!r} to type."
                    ).format(opt, str)
                )
    else:
        raise \
            optparse.OptionValueError(
                "option {0!s}: invalid list value: {1!r}".format(opt, str)
            )

def checkIntList(option, opt, value):
    """
    Extracts a list of integers from a specified value.
    @param option: The option.
    @type opt: string
    @param opt: The command line option string.
    @type value: string
    @param value: The value passed to the option.
    @rtype list
    @return: list of ints
    """
    return getListFromString(opt, value, int)

def checkFloatList(option, opt, value):
    """
    Extracts a list of floats from a specified value.
    @type opt: string
    @param opt: The command line option string.
    @type value: string
    @param value: The value passed to the option.
    @rtype list
    @return: list of floats
    """
    return getListFromString(opt, value, float)

class ImprovedOption(optparse.Option):
    """
    Extends optparse.Option class by adding 'int_list' and 'float_list' types
    """
    TYPES = optparse.Option.TYPES + ("int_list", "float_list",)
    TYPE_CHECKER = copy.copy(optparse.Option.TYPE_CHECKER)
    TYPE_CHECKER["int_list"] = checkIntList
    TYPE_CHECKER["float_list"] = checkFloatList

class OptionParser(optparse.OptionParser):
    """
    Command line option parser for use in conjunction with
    simple wave propagation simulations.
    """
    def __init__(self, usage=None):
        """
        Initialises this parser with options.
        @type usage: string
        @param usage: Usage string, if None, default value is given.
        """
        if (usage==None):
            usage =\
              "usage: %prog [options]\n\n" +\
              "Runs wave propagation simulation. A rectangular block of" +\
              "\nparticles acts as the elastic medium. Two different spring" +\
              "\nconstants can be specified to create non-homogeneity." +\
              "\nA single point source disturbance creates a wave in the" +\
              "elastic block."
        optparse.OptionParser.__init__(
            self,
            usage=usage,
            option_class=ImprovedOption
        )
        
        self.add_option(
          "-n", "--mpi-worker-processes",
          dest="numWorkerProcesses",
          type="int",
          metavar="N",
          default=2,
          help=\
              "Defines the number of MPI worker processes " +\
              "(default N=2)"
        )

        self.add_option(
          "-d", "--mpi-dim-list",
          dest="mpiDimList",
          type="int_list",
          metavar="D",
          default=[2,0,0],
          help=\
              "Defines the spatial grid division of domain among "+\
              " MPI worker processes (default D=[2,0,0])"
        )

        self.add_option(
          "-t", "--time-step-size",
          dest="timeStepSize",
          type="float",
          metavar="T",
          default=0.05,
          help=\
              "Defines time step size used in explicit integration scheme "+\
              "(default T=0.05)"
        )

        self.add_option(
          "-r", "--radius",
          dest="radius",
          type="float",
          metavar="R",
          default=1.0,
          help=\
              "Defines radius of particles "+\
              "(default R=1.0)"
        )

        self.add_option(
          "-m", "--max-time-steps",
          dest="maxNumTimeSteps",
          type="int",
          metavar="M",
          default=2000,
          help=\
              "Defines max number of time step iterations "+\
              "(default M=2000)"
        )

        self.add_option(
          "-b", "--block-size",
          dest="particlesPerDim",
          type="int_list",
          metavar="B",
          default=[128,128,1],
          help=\
              "Defines the size of the particle block as number of " +\
              "particles per dimension component "+\
              "(default B=[128,128,1])"
        )

        self.add_option(
          "-P", "--particle-data-incr",
          dest="particleDataIncr",
          type="int",
          metavar="N",
          default=100,
          help=\
              "Defines the frequency at which particle displacement data is " +\
              "saved to file. Data is saved every N time steps "+\
              "(default N=100)"
        )

        self.add_option(
          "-S", "--seismo-data-incr",
          dest="seismoDataIncr",
          type="int",
          metavar="N",
          default=10,
          help=\
              "Defines the frequency at which seismograph data is " +\
              "saved to file. Data is saved every N time steps "+\
              "(default N=10)"
        )

        self.add_option(
          "-v", "--verbosity",
          dest="verbosity",
          action="store_true",
          default=False,
          help=\
              "Generates (lots of) debug output during simulation run" +\
              " (default no debug output)"
        )

        self.add_option(
          "-s", "--source-depth",
          dest="sourceDepth",
          type="float",
          default=0.5,
          metavar="D",
          help=\
              "Specifies the depth of the source disturbance as a value" +\
              " between 0.0 and 1.0. A value of 0.0 places the source at" +\
              " the surface of the particle-block, a value of 0.5 places the" +\
              " source at the centre of the block, a value of 1.0 places" +\
              " places the source at the bottom of the particle-block, etc"
              " (default D=0.5). The source is centred in the particle block" +\
              " with respect to the x and z dimensions."
        )

        self.add_option(
          "-f", "--source-frequency",
          dest="sourceFrequency",
          type="float",
          default=0.02,
          metavar="F",
          help=\
              "Specifies the frequency of the sinusoidal source disturbance" +\
              " (default F=0.02)"
        )

        self.add_option(
          "-D", "--source-max-displacement",
          dest="sourceMaxDisplacement",
          type="float_list",
          default=[0.1,0.1,0.0],
          metavar="M",
          help=\
              "Specifies the maximum relative displacement point of" +\
              " the source disturbance "+\
              " (default F=[0.1,0.1,0.0]). If the point-source is located at" +\
              " point X, then the source particle will be moved on the line" +\
              " between the point X and the point X+F."
        )

        self.add_option(
          "-U", "--upper-medum-depth",
          dest="upperMediumDepth",
          type="float",
          metavar="D",
          default=0.25,
          help=\
              "Specifies the depth of the --upper-spring-K bonds as a value" +\
              " between 0.0 and 1.0. A 0.0 value implies no upper-spring-K" +\
              " bonds, a 0.5 value creates half the bonds as upper-spring-K" +\
              " and half the bonds as lower-spring-K, a 1.0 value "
              " creates all bonds as upper-spring-K, etc" +\
              " (default D=0.25)"
        )

        self.add_option(
          "-u", "--upper-spring-K",
          dest="upperSpringK",
          type="float",
          metavar="K",
          default=1.0,
          help=\
              "Defines bond spring constant for 'upper' particle region"+\
              " (default K=1.0)"
        )
        self.add_option(
          "-l", "--lower-spring-K",
          dest="lowerSpringK",
          type="float",
          metavar="K",
          default=1.0,
          help=\
              "Defines bond spring constant for 'lower' particle region"+\
              " (default K=1.0)"
        )
