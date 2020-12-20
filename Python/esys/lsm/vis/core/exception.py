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
Defines helper function L{raiseNotImplemented} for raising an exception
in abstract methods of base classes.
"""
import esys.lsm.Logging
import logging
import traceback
import sys

_defaultLogger = esys.lsm.Logging.getLogger("esys.lsm.vis.core")

_doRaiseWhenNotImplemented = True

def raiseExceptionWhenNotImplemented(doRaise=None):
    """
    Return value of this function used in L{raiseNotImplemented}.
    @type doRaise: bool
    @param doRaise: If specified, sets whether L{raiseNotImplemented}
    is to raise an exception.
    @rtype: bool
    @return: True when L{raiseNotImplemented} is to raise
    a C{NotImplementedError} exception.
    """
    global _doRaiseWhenNotImplemented
    if (doRaise != None):
        _doRaiseWhenNotImplemented = doRaise
    return _doRaiseWhenNotImplemented

def raiseNotImplemented(msg = "", logger = _defaultLogger):
    """
    Raises a NotImplementedError exception and logs
    the back-trace using the specified logger. The exception
    is raised only if L{raiseExceptionWhenNotImplemented}
    returns true, otherwise only the back-trace is logged
    and no exception is thrown.
    @type msg: str
    @param msg: Exception message string passed to NotImplementedError
    constructor.
    @type logger: C{logging.Logger}
    @param logger: The C{Logger.error} message is used to log
    the lines of a back-trace.
    """
    try:
        raise NotImplementedError(msg)
    except:
        if (logger != None):
            traceBackMsg = \
                traceback.format_exception(
                    sys.exc_info()[0],
                    sys.exc_info()[1],
                    sys.exc_info()[2]
                )
            for line in traceBackMsg:
                for subLine in str.split(str.strip(line), "\n"):
                    logger.error(subLine)

        if (raiseExceptionWhenNotImplemented()):
            raise
