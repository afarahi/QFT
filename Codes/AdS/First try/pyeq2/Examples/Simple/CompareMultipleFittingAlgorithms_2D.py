#    Version info: $Id: CompareMultipleFittingAlgorithms_2D.py 28 2012-05-12 13:13:18Z zunzun.com@gmail.com $

import os, sys, inspect

# ensure pyeq2 can be imported
if os.path.join(sys.path[0][:sys.path[0].rfind(os.sep)], '../..') not in sys.path:
    sys.path.append(os.path.join(sys.path[0][:sys.path[0].rfind(os.sep)], '../..'))
import pyeq2



if __name__ == "__main__":
    fittingTargetText = 'SSQABS'
    
    deEstimatedCoefficients = []
    
    print 'It is very rare for an algorithm to fit better than Levenberg-Marquardt,'
    print 'however this particular example shows that it is at least barely possible.'
    print
    
    for fittingAlgorithmName in pyeq2.solverService.ListOfNonLinearSolverAlgorithmNames[:2]: # only two algorithms used in this example
        equation = pyeq2.Models_2D.BioScience.AphidPopulationGrowth(fittingTargetText, 'Offset')
        

        if equation.CanLinearSolverBeUsedForSSQABS() == True and fittingTargetText == 'SSQABS':
            raise Exception('The selected combination of equation and SSQABS fitting target does not use a non-linear solver')
        
        if fittingTargetText == 'ODR':
            raise Exception('ODR cannot use multiple fitting algorithms')        

        
        pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(equation.exampleData, equation, False)
        
        equation.deEstimatedCoefficients = deEstimatedCoefficients
        equation.Solve(inNonLinearSolverAlgorithmName=fittingAlgorithmName)
        deEstimatedCoefficients = equation.deEstimatedCoefficients # no need to re-run DE
        
        print fittingTargetText, 'of', equation.CalculateAllDataFittingTarget(equation.solvedCoefficients), 'for the fitting algorithm', fittingAlgorithmName
        print 'Coefficients:', equation.solvedCoefficients
        print
        sys.stdout.flush()
