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
Defines L{OptionParser} class which extends C{optparse.OptionParser} class
by adding 'string_list', 'int_list' and 'float_list' option types. Also
defines L{LogOptionParser} which extends L{OptionParser} and adds support
for L{esys.lsm.Logging} options.
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
            return list(map(mapCallable, str.split(",", match.group(1))))
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
    @rtype: list
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
    @rtype: list
    @return: list of floats
    """
    return getListFromString(opt, value, float)

def checkStringList(option, opt, value):
    """
    Extracts a list of strings from a specified value.
    @type opt: string
    @param opt: The command line option string.
    @type value: string
    @param value: The value passed to the option.
    @rtype: list
    @return: list of floats
    """
    return getListFromString(opt, value, str)

class ListOption(optparse.Option):
    """
    Extends optparse.Option class by adding 'string_list',
    'int_list' and 'float_list' types
    """
    TYPES = optparse.Option.TYPES + ("int_list", "float_list", "string_list")
    TYPE_CHECKER = copy.copy(optparse.Option.TYPE_CHECKER)
    TYPE_CHECKER["int_list"] = checkIntList
    TYPE_CHECKER["float_list"] = checkFloatList
    TYPE_CHECKER["string_list"] = checkStringList

class OptionParser(optparse.OptionParser):
    """
    Command line option parser which extends C{optparse.OptionParser}
    by adding "int_list", "float_list" and "string_list" types.
    """
    def __init__(
        self,
        usage=None,
        option_list=None,
        option_class=ListOption,
        version=None,
        conflict_handler="error",
        description=None,
        formatter=None,
        add_help_option=1,
        prog=None
    ):
        """
        Initialises this parser, arguments as per C{optparse.OptionParser}.
        """
        optparse.OptionParser.__init__(
            self,
            usage=usage,
            option_list=option_list,
            option_class=option_class,
            version=version,
            conflict_handler=conflict_handler,
            description=description,
            formatter=formatter,
            add_help_option=add_help_option,
            prog=prog
        )

from esys.lsm import Logging

class LogOptionParser(OptionParser):
    """
    Command line option parser which extends L{OptionParser}
    by adding a default logging option.
    """
    def __init__(
        self,
        usage=None,
        option_list=None,
        option_class=ListOption,
        version=None,
        conflict_handler="error",
        description=None,
        formatter=None,
        add_help_option=1,
        prog=None
    ):
        """
        Initialises this parser, arguments as per L{OptionParser}.
        """
        OptionParser.__init__(
            self,
            usage=usage,
            option_list=option_list,
            option_class=option_class,
            version=version,
            conflict_handler=conflict_handler,
            description=description,
            formatter=formatter,
            add_help_option=add_help_option,
            prog=prog
        )
        self.addDefaultOptions()

    def addDefaultOptions(self):
        self.add_option(
          "-l", "--log-level",
          dest="logLevel",
          type="string",
          default="WARNING",
          metavar="L",
          help=\
              "The level of logging output (default L=\"WARNING\")"
        )

    def initialiseLogging(self, options):
        logLevel = options.logLevel
        if ((logLevel != None) and (len(logLevel) > 0)):
            Logging.basicConfig()
            Logging.getLogger("").setLevel(Logging.getLevel(logLevel))

    def parse_args (self, args=None, values=None):
        (options,values) = OptionParser.parse_args(self, args, values)
        self.initialiseLogging(options)
        return (options,values)
