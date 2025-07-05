"""A collection of latex printing code."""

from .utility import sumColumns
from functools import cmp_to_key

def beginDoc():
    return "\\documentclass{article}\n" + \
        "\\usepackage{tikz}\n" + \
        "\\begin{document}\n"

def endDoc():
    return "\\end{document}\n"

def beginTikz(scale):
    return "\\begin{tikzpicture} [x=%dpt, y=%dpt," % (scale, scale) + \
        "baseline=(current bounding box.center)]\n"

def endTikz():
    return "\\end{tikzpicture}\n"

def showLocalStrandDiagram(sd, pmc_map = None):
    """sd is the local strand diagram to be displayed. pmc_map is a mapping from
    points on the local PMC to vertical location in the displayed picture.
    Defaults to i -> i for all 0 <= i < n.

    """
    result = ""
    local_pmc = sd.local_pmc
    n = local_pmc.n
    if pmc_map is not None:
        assert len(pmc_map) == n
    else:
        pmc_map = [i for i in range(n)]

    # Header
    result += beginTikz(10)

    # Display a dummy line, to align the strand diagrams and idempotents
    min_pos, max_pos = min(pmc_map), max(pmc_map)
    result += "\\draw [color=white] (3, %d) to (3, %d);\n" \
              % (min_pos - 1, max_pos + 1)

    # Display dots for the local PMC.
    for i in range(n):
        if i not in local_pmc.endpoints:
            result += "\\filldraw (0, %d) circle (0.1);\n" % pmc_map[i]
            result += "\\filldraw (3, %d) circle (0.1);\n" % pmc_map[i]

    # Display single and double horizontal lines.
    for idem in sd.all_hor:
        for p in local_pmc.pairs[idem]:
            result += "\draw [dashed] (0, %d) to (3, %d);\n" \
                      % (pmc_map[p], pmc_map[p])

    # Display strands.
    for start, end in sd.strands:
        result += "\draw [->] (0, %d) to [out=0, in=180] (3, %d);\n" \
                  % (pmc_map[start], pmc_map[end])

    # Footer
    result += endTikz()

    return result

def showDAGenerator(gen, pmc_map1, pmc_map2):
    """Show a generator of the type DA structure."""
    result = ""

    # Header
    result += beginTikz(10)

    local_pmc1 = gen.idem1.local_pmc
    local_pmc2 = gen.idem2.local_pmc
    n1 = local_pmc1.n
    n2 = local_pmc2.n

    # Display local PMC dots and middle vertible line.
    for i in range(n1):
        if i not in local_pmc1.endpoints:
            result += "\\filldraw (0, %d) circle (0.1);\n" % pmc_map1[i]
    for i in range(n2):
        if i not in local_pmc2.endpoints:
            result += "\\filldraw (6, %d) circle (0.1);\n" % pmc_map2[i]
    min_pos = min(min(pmc_map1), min(pmc_map2))
    max_pos = max(max(pmc_map1), max(pmc_map2))
    result += "\\draw (3, %d) to (3, %d);\n" % (min_pos - 1, max_pos + 1)

    # Display single and double horizontals
    for idem in gen.idem1:
        for p in local_pmc1.pairs[idem]:
            result += "\\draw [dashed] (0, %d) to (3, %d);\n" \
                      % (pmc_map1[p], pmc_map1[p])
    for idem in gen.idem2:
        for p in local_pmc2.pairs[idem]:
            result += "\\draw [dashed] (3, %d) to (6, %d);\n" \
                      % (pmc_map2[p], pmc_map2[p])

    # Footer
    result += endTikz()
    return result

def showDAStructure(dastr, pmc_map1 = None, pmc_map2 = None):
    """Show a type DA structure."""
    result = ""

    n1 = dastr.algebra1.local_pmc.n
    n2 = dastr.algebra2.local_pmc.n
    if pmc_map1 is not None:
        assert len(pmc_map1) == n1
    else:
        pmc_map1 = [i for i in range(n1)]
    if pmc_map2 is not None:
        assert len(pmc_map2) == n2
    else:
        pmc_map2 = [i for i in range(n2)]

    to_print = []

    for (gen_from, coeffs_a), target in list(dastr.da_action.items()):
        for (coeff_d, gen_to), ring_coeff in list(target.items()):
            to_print.append((gen_from, coeffs_a, coeff_d, gen_to))

    # Sort the to_print by the following comparison function
    def compare_arrow(arrow1, arrow2):
        gen_from1, coeffs_a1, coeff_d1, gen_to1 = arrow1
        gen_from2, coeffs_a2, coeff_d2, gen_to2 = arrow2
        def compare_mult(mult1, mult2):
            if sum(mult1) != sum(mult2):
                return sum(mult1) - sum(mult2)
            mult1, mult2 = list(reversed(mult1)), list(reversed(mult2))
            if mult1 < mult2:
                return -1
            if mult1 > mult2:
                return 1
            return 0

        mult_d1, mult_d2 = coeff_d1.multiplicity, coeff_d2.multiplicity
        if compare_mult(mult_d1, mult_d2) != 0:
            return compare_mult(mult_d1, mult_d2)
        if len(coeff_d1.strands) != len(coeff_d2.strands):
            return len(coeff_d1.strands) - len(coeff_d2.strands)

        mult_a1 = sumColumns([coeff.multiplicity for coeff in coeffs_a1], n1-1)
        mult_a2 = sumColumns([coeff.multiplicity for coeff in coeffs_a2], n2-1)
        return compare_mult(mult_a1, mult_a2)

    for gen_from, coeffs_a, coeff_d, gen_to in sorted(to_print,
                                                      key = cmp_to_key(compare_arrow)):
        result += "\\begin{equation}\n"
        result += "\\delta^1\\left( \n"
        result += showDAGenerator(gen_from, pmc_map1, pmc_map2)
        for coeff_a in coeffs_a:
            result += "~,~\n";
            result += showLocalStrandDiagram(coeff_a, pmc_map2)
        result += "\\right) \\to \n"
        result += showLocalStrandDiagram(coeff_d, pmc_map1)
        result += "~\\otimes~ \n"
        result += showDAGenerator(gen_to, pmc_map1, pmc_map2)
        result += "\\end{equation}\n"

    return result

def showArrow(coeff_d, coeffs_a, pmc_map1 = None, pmc_map2 = None):
    result = ""
    result += "\\begin{equation}\n"
    result += "~,~\n".join([showLocalStrandDiagram(coeff_a, pmc_map2)
                            for coeff_a in coeffs_a])
    result += "\\to \n"
    result += showLocalStrandDiagram(coeff_d, pmc_map1)
    result += "\\end{equation}\n"

    return result
