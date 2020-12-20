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
Defines classes which relate to I{modifying} renderable objects.
"""

class Modifiable(object):
    def __init__(self):
        self.idModDict = dict()

    def getModifierList(self):
        pairList = list(self.idModDict.values())
        pairList.sort(key=lambda a: a[0])
        return [e[1] for e in pairList]

    def applySingle(self, modifier):
        if (modifier != None):
            self.idModDict[id(modifier.__class__)] = \
              (len(self.idModDict),modifier)

    def applySequence(self, sequence):
        for modifier in sequence:
            self.applySingle(modifier)

    def apply(self, modifiers):
        if ((modifiers != None) and (hasattr(modifiers, "__iter__"))):
            self.applySequence(modifiers)
        else:
            self.applySingle(modifiers)

class Modifier(object):
    def __init__(self):
        pass

