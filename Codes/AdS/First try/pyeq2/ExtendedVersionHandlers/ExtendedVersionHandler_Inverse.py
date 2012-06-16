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
#    Version info: $Id: ExtendedVersionHandler_Inverse.py 14 2012-02-16 11:13:54Z zunzun.com@gmail.com $

import pyeq2
import IExtendedVersionHandler


class ExtendedVersionHandler_Inverse(IExtendedVersionHandler.IExtendedVersionHandler):
    
    def AssembleDisplayHTML(self, inModel):
        if inModel.GetDimensionality() == 2:
            return inModel._HTML + '<br>' + inModel._leftSideHTML + ' = x / ' + inModel._leftSideHTML
        else:
            return inModel._HTML + '<br>' + inModel._leftSideHTML + ' = xy / ' + inModel._leftSideHTML


    def AssembleDisplayName(self, inModel):
        return 'Inverse ' + inModel._baseName


    def AssembleSourceCodeName(self, inModel):
        return inModel.__class__.__name__ + "_Inverse"


    def AssembleCoefficientDesignators(self, inModel):
        return inModel._coefficientDesignators


    # overridden from abstract parent class
    def AppendAdditionalCoefficientBounds(self, inModel):
        return


    def AssembleOutputSourceCodeCPP(self, inModel):
        if inModel.GetDimensionality() == 2:
            return inModel.SpecificCodeCPP() + "\ttemp = x_in / temp;\n"
        else:
            return inModel.SpecificCodeCPP() + "\ttemp = (x_in * y_in) / temp;\n"


    def GetAdditionalModelPredictions(self, inBaseModelCalculation, inCoeffs, inDataCacheDictionary, inModel):
        if inModel.GetDimensionality() == 2:
            return self.ConvertInfAndNanToLargeNumber(inDataCacheDictionary['X'] / inBaseModelCalculation)
        else:
            return self.ConvertInfAndNanToLargeNumber(inDataCacheDictionary['XY'] / inBaseModelCalculation)
