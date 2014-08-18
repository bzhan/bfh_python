"""Producing type DA structures for Dehn twists, using local actions."""

from autocompleteda import AutoCompleteDAStructure
from dastructure import DAStructure, SimpleDAGenerator
from extendbyid import ExtendedDAStructure, LocalDAStructure
from localpmc import LocalIdempotent, LocalStrandAlgebra, PMCSplitting
from pmc import linearPMC
from utility import subset
from utility import F2
import ast
import itertools

class AntiBraidDA(ExtendedDAStructure):
    """Responsible for producing a type DA structure for the anti-braid
    resolution (admissible case only), using local actions.

    """
    def __init__(self, genus, c_pair):
        """Specifies genus of the starting PMC and the ID of the pair of
        anti-braid resolution.

        """
        self.genus = genus
        self.c_pair = c_pair
        self.pmc = linearPMC(genus)

        self.n = 4 * genus
        self.c1, self.c2 = self.pmc.pairs[c_pair]

        if self.c2 == self.c1 + 3:
            self.is_degenerate = False
        else:
            assert self.c2 == self.c1 + 2
            assert self.c1 == 0 or self.c2 == self.n - 1
            self.is_degenerate = True

        if self.is_degenerate:
            # One position between c1 and c2, called p
            self.p = self.c1 + 1
            self.p_pair = self.pmc.pairid[self.p]
        else:
            # Two positions between c1 and c2, for (d)own and (u)p
            self.d = self.c1 + 1
            self.u = self.c1 + 2
            self.d_pair = self.pmc.pairid[self.d]
            self.u_pair = self.pmc.pairid[self.u]

        # Necessary to get local DA structure.
        self.splitting = PMCSplitting(self.pmc, [(self.c1, self.c2)])
        self.local_pmc = self.splitting.local_pmc
        self.mapping = self.splitting.local_mapping

        # Local DA Structure
        self.local_da = self.getLocalDAStructure()

        ### Uncomment to use autocompleteda to construct arrows from seeds.
        # auto = AutoCompleteDAStructure(self.local_da, ([]))
        # self.local_da = auto.complete()

        # Initiate the ExtendedDAStructure
        ExtendedDAStructure.__init__(self, self.local_da,
                                     self.splitting, self.splitting)

    def getLocalDAStructure(self):
        """Returns the local type DA structure associated to the anti-braid
        resolution.

        """
        if self.is_degenerate:
            if self.c1 == 0:
                patterns_raw = self._get_patterns_bottom()
            else:
                assert self.c2 == self.n - 1
                patterns_raw = self._get_patterns_top()
        else:
            patterns_raw = self._get_patterns_middle()

        arrow_patterns = {}
        for pattern in patterns_raw:
            start_class, end_class = pattern[0], pattern[1]
            coeffs_a = []
            for i in range(2, len(pattern)-1):
                coeffs_a.append(self.local_pmc.sd(pattern[i]))
            key = (start_class, end_class, tuple(coeffs_a))
            if key not in arrow_patterns:
                arrow_patterns[key] = []
            arrow_patterns[key].append(self.local_pmc.sd(pattern[-1]))

        # Now start construction of the local DA structure.
        alg = LocalStrandAlgebra(F2, self.local_pmc)

        local_c_pair = 0
        # Compute the set of local generators. Generators of class 0 has
        # idempotents (l_idem, r_idem) where l_idem has the c_pair and
        # r_idem has one of the u or d pairs (with the rest being the same).
        # Generators of class 1 and 2 has l_idem = r_idem such that c_pair
        # is in both.
        if self.is_degenerate:
            single_idems = [1]  # local p_pair
            da_idems_0 = [([0], [1])]
        else:  # Non-degenerate case
            single_idems = [1, 2]  # local d_pair and local u_pair
            da_idems_0 = [([0], [1]), ([0], [2]),
                          ([0, 1], [1, 2]), ([0, 2], [1, 2])]  # class 0

        local_da = LocalDAStructure(
            F2, alg, alg, single_idems1 = single_idems,
            single_idems2 = single_idems)
        for i in range(len(da_idems_0)):  # class 0
            l_idem, r_idem = da_idems_0[i]
            local_da.addGenerator(SimpleDAGenerator(
                local_da, LocalIdempotent(self.local_pmc, l_idem),
                LocalIdempotent(self.local_pmc, r_idem), "0_%d" % i))
        all_idems = subset(range(self.local_pmc.num_pair))  # class 1 and 2
        for i in range(len(all_idems)):
            idem = LocalIdempotent(self.local_pmc, all_idems[i])
            if local_c_pair in idem:
                local_da.addGenerator(
                    SimpleDAGenerator(local_da, idem, idem, "1_%d" % i))
                local_da.addGenerator(
                    SimpleDAGenerator(local_da, idem, idem, "2_%d" % i))

        mod_gens = local_da.getGenerators()
        # Have to take care of u_map. It is sufficient to know that u_map musst
        # preserve class of generators.
        for i in range(len(single_idems)):
            idem = single_idems[i]
            for local_gen in mod_gens:
                idem1, idem2 = local_gen.idem1, local_gen.idem2
                if idem in idem1 and idem in idem2:
                    # local_gen is eligible for u_maps[i]
                    target_idem1 = idem1.removeSingleHor([idem])
                    target_idem2 = idem2.removeSingleHor([idem])
                    target_gen = [target for target in mod_gens
                                  if target.idem1 == target_idem1 and
                                  target.idem2 == target_idem2 and
                                  target.name[0] == local_gen.name[0]]
                    assert len(target_gen) == 1
                    local_da.add_u_map(i, local_gen, target_gen[0])
        # Check all u_map has been filled.
        local_da.auto_u_map()

        # Add arrows according to arrow_pattern.
        for key in arrow_patterns.keys():
            start_class, end_class, coeffs_a = key
            if len(coeffs_a) == 1 and coeffs_a[0].isIdempotent():
                continue
            for coeff_d in arrow_patterns[key]:
                used = False
                for x, y in itertools.product(mod_gens, mod_gens):
                    if x.name[0] == "%d" % start_class and \
                       y.name[0] == "%d" % end_class and \
                       DAStructure.idemMatchDA(x, y, coeff_d, coeffs_a):
                        local_da.addDelta(x, y, coeff_d, coeffs_a, 1)
                        used = True
                if not used:
                    print "Warning: unused arrow: %s %s" % (coeffs_a, coeff_d)
        return local_da

    def _get_patterns_middle(self):
        """Returns the local patterns."""
        # Local PMC is 0*-1-2-3-4-5*, with 1 and 4 paired.
        input_patterns = open("antibraid_arrows.data", "r")
        patterns_raw = ast.literal_eval(input_patterns.read())
        return patterns_raw

    def _get_patterns_bottom(self):
        # Local PMC is 0-1-2-3*, with 0 and 2 paired.
        patterns_raw = [
            # Initial patterns
            (1, 2, [0]), (1, 2, [0, 1]),

            (0, 2, [(1, 2)], [0]),
            (1, 0, [(0, 1)], [0]),
            (1, 2, [(0, 2)], [0]),
            (1, 2, [1, (0, 2)], [0, 1]),

            # Seed for the middle regions
            (0, 0, [(0, 2)]),

            # Added for the middle regions
            (2, 2, [(0, 2)]),
            (2, 2, [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
            (1, 1, [(0, 2)]),
            (1, 1, [1, (0, 2)]),
            (1, 2, [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
            (2, 1, [(0, 1), (1, 2)]),
            (1, 1, [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
            (2, 2, [1, (0, 2)]),
        ]
        return patterns_raw

    def _get_patterns_top(self):
        # Local PMC is 0*-1-2-3, with 1 and 3 paired.
        # Simply add one to everything in _get_patterns_bottom
        def translate(pattern):
            result = []
            for entry in pattern:
                if isinstance(entry, int):
                    result.append(entry + 1)
                else:
                    result.append((entry[0] + 1, entry[1] + 1))
            return result

        patterns_raw = []
        for arrow_pattern in self._get_patterns_bottom():
            patterns_raw.append(
                arrow_pattern[:2] + tuple([translate(pattern)
                                           for pattern in arrow_pattern[2:]]))
        return patterns_raw
