#    pyeq2 is a collection of equations expressed as Python classes
#
#    Copyright (C) 2012 James R. Phillips
#    2548 Vera Cruz Drive
#    Birmingham, AL 35235 USA
#
#    email: zunzun@zunzun.com
#    web: http://zunzun.com
#
#    License: BSD-style (see LICENSE.txt in main source directory)
#    Version info: $Id: ExtendedVersionHandler_Offset.py 21 2012-03-10 19:48:51Z zunzun.com@gmail.com $

import pyeq2
import IExtendedVersionHandler


class ExtendedVersionHandler_Offset(IExtendedVersionHandler.IExtendedVersionHandler):
    
    def AssembleDisplayHTML(self, inModel):
        return inModel._HTML + " + Offset"


    def AssembleDisplayName(self, inModel):
        return inModel._baseName + " With Offset"


    def AssembleSourceCodeName(self, inModel):
        return inModel.__class__.__name__ + "_Offset"


    def AssembleCoefficientDesignators(self, inModel):
        return inModel._coefficientDesignators + ['Offset']


    # overridden from abstract parent class
    def AppendAdditionalCoefficientBounds(self, inModel):
        if inModel.upperCoefficientBounds != []:
            inModel.upperCoefficientBounds.append(None)
        if inModel.lowerCoefficientBounds != []:
            inModel.lowerCoefficientBounds.append(None)


    def AssembleOutputSourceCodeCPP(self, inModel):
        return inModel.SpecificCodeCPP() + "\ttemp += Offset;\n"


    # overridden from abstract parent class
    def GetAdditionalDataCacheFunctions(self, inModel, inDataCacheFunctions):
        return inDataCacheFunctions


    def GetAdditionalModelPredictions(self, inBaseModelCalculation, inCoeffs, inDataCacheDictionary, inModel):
        return self.ConvertInfAndNanToLargeNumber(inBaseModelCalculation + inCoeffs[len(inCoeffs)-1])
    
    
    # overridden from abstract parent class
    def CanLinearSolverBeUsedForSSQABS(self, inModelFlag):
        return False

