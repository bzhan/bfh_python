"""Sign conventions."""

from fractions import Fraction
from grading import standardRefinement
from pmc import StrandDiagram

class AbsZ2Grading:
    """Computes absolute Z/2Z grading for the strand algebra."""
    def __init__(self, pmc):
        """Creates data needed for computation of Z/2Z grading for the given
        pointed matched circle.

        """
        self.pmc = pmc

    def getAbsGrading(self, sd):
        """Returns the absolute Z/2Z grading (returned value is either 0 or 1).

        """
        small_gr = sd.getSmallGrading(refinement = standardRefinement);
        gr_group = small_gr.parent
        for i in range(len(small_gr.spinc)):
            mult = -small_gr.spinc[i]
            small_gr = small_gr * gr_group.basis(i).power(mult)
            small_gr.maslov += Fraction(-1, 2) * mult
        assert small_gr.maslov.denominator == 1
        return small_gr.maslov.numerator % 2
