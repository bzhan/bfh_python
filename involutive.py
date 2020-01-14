"""
Created on Wed May 31 11:22:54 2017

@author: lipshitz
"""
from utility import SummableDict, F2, fracToInt, ACTION_LEFT, ACTION_RIGHT
from algebra import Element, E0
from dstructure import MorDtoDGenerator, DGenerator, SimpleDStructure
from algebra import SimpleChainComplex, SimpleGenerator, SimpleChainMorphism, Generator
from dastructure import SimpleDAStructure, SimpleDAGenerator, identityDA, DATensorDGenerator
from pmc import StrandDiagram
from grading import SmallGradingGroup, SimpleDbGradingSet, SimpleDbGradingSetElement
from braid import Braid, BraidCap, readBridgePresentation
from arcslideda import ArcslideDA
    

def applyDMor(f,x):
    "Apply a type D morphism, of Element class, to a DGenerator x"
    answer = E0
    for f0 in list(f.keys()):
        answer += f[f0]*f0.apply(x)
    return answer

def composeMor(f,g, parent=None):
    "Return the composition f\circ g, if f and g are type D morphism Elements (gotten, say, from prev_meaning)"
    answer = E0
    for g0 in list(g.keys()): #g0 is of type WHAT?
        for f0 in list(f.keys()):
            addit = (f[f0]*g[g0])* f0.compose(g0,parent)
            answer += addit
    return answer


def chordPairs(pmc):
    "List of chord pairs as in differential on CFDD(Id)."
    answer = list()
    algebra = pmc.getAlgebra()
    for i in range(pmc.n):
        for j in range(i+1,pmc.n):
            for I in pmc.getIdempotents():
                Ipairs = [I.pmc.pairs[x] for x in I]
                occupied = list()
                for x in Ipairs:
                    occupied.extend(x) #occupied positions in I
                if (i in occupied) and not (j in occupied):
                    sigma = StrandDiagram(algebra,I,[(i,j)])
                    J = sigma.right_idem.comp()
                    rho = StrandDiagram(algebra,J,[(i,j)])
                    answer.append((sigma,rho))
    return answer

def azDA(pmc,mult_one = True):
    "The type DA module associated to the Auroux-Zarev piece. Input: the alpha pointed matched circle."
    algebra = pmc.getAlgebra(mult_one = mult_one)
    answer = SimpleDAStructure(F2, algebra, algebra)
    #Compute the generators
    for a in pmc.getStrandDiagrams(answer.algebra1):
        gen = SimpleDAGenerator(answer, a.left_idem.comp(), a.right_idem, a)
        answer.addGenerator(gen)
    #Terms in the differential coming from the differential on the algebra:
    for x in answer.generators:
        dx = x.name.diff()
        for y in list(dx.keys()):
            ygen = SimpleDAGenerator(answer, y.left_idem.comp(), y.right_idem, y)
            answer.addDelta(x,ygen,x.idem1.toAlgElt(algebra),list(),dx[y])
    #Terms coming from multiplying by algebra elements on the right
    for x in answer.generators:
        xa = x.name
        could_multiply = algebra.getGeneratorsForIdem(left_idem = xa.right_idem)
        for b in could_multiply:
            if not b.isIdempotent():
                xab = xa*b
                for y in list(xab.keys()):
                    ygen = SimpleDAGenerator(answer,y.left_idem.comp(), y.right_idem, y)
                    answer.addDelta(x,ygen,x.idem1.toAlgElt(algebra),[b,],xab[y])
    #Terms coming from differential on CFDD(Id):
    chord_pairs = chordPairs(pmc)
    for (sigma, rho) in chord_pairs:
        for x in algebra.getGeneratorsForIdem(left_idem = rho.right_idem):
            xgen = SimpleDAGenerator(answer,x.left_idem.comp(), x.right_idem, x)
            rhox = rho*x
            for y in list(rhox.keys()):
                ygen = SimpleDAGenerator(answer, y.left_idem.comp(), y.right_idem, y)
                answer.addDelta(xgen, ygen, sigma, list(),rhox[y])    
    return answer


def tensorDAid(P,Q,M,MP,MQ,morcx,f):
    """Input: 
        --simple type D structures P and Q
        --a simple type DA structure M
        --the tensor products M\boxtimes P and M\boxtimes Q
        --A morphism f in Mor(P,Q).
    Returns: Id\boxtimes f, where Id is the identity morphism of M.
    """
    answer = E0    

    def search(start_gen, cur_dgen, cur_coeffs_a, inTarget = False):
        """Searching for an arrow in the box tensor product.
        - start_gen: starting generator in the box tensor product. The
        resulting arrow will start from here.
        - cur_dgen: current location in the type D structure.
        - cur_coeffs_a: current list of A-side inputs to the type DA
        structure (or alternatively, list of algebra outputs produced by
        the existing path through the type D structure).
        """
        answer = E0
        start_dagen, start_dgen = start_gen
        cur_delta = M.delta(start_dagen, cur_coeffs_a)
        for (coeff_d, gen_to), ring_coeff in list(cur_delta.items()):
            if inTarget:
                startelt = DATensorDGenerator(MP,start_dagen,start_dgen)
                endelt = DATensorDGenerator(MQ,gen_to,cur_dgen)
                morph = MorDtoDGenerator(morcx, startelt, coeff_d, endelt)
                assert morph.source.parent == MP and morph.target.parent == MQ
                answer += 1*morph
        if M.deltaPrefixNS(start_dagen, cur_coeffs_a):
            if not inTarget:
                for (coeff_out, dgen_to), ring_coeff in list(P.delta(cur_dgen).items()):
                    answer += search(start_gen, dgen_to, cur_coeffs_a + (coeff_out,), inTarget=False)
                for (coeff_out, dgen_to), ring_coeff in list(applyDMor(f,cur_dgen).items()): #APPLY f
                    answer += search(start_gen, dgen_to, cur_coeffs_a + (coeff_out,), inTarget=True)
            if inTarget:
                for (coeff_out, dgen_to), ring_coeff in list(Q.delta(cur_dgen).items()):
                    answer += search(start_gen, dgen_to, cur_coeffs_a + (coeff_out,), inTarget=True)
        return answer

    for x in MP.getGenerators():
        #Invoke search here with starting data coming from x.
        dagen, dgen = x
        answer += search(x, dgen, ())

    return answer

def mappingConeD(f, returnDicts = False):
    "Return the mapping cone of a morphism f (Element class) from one SimpleDStructure to another."
    #Extract some basic info from f
    for f0 in list(f.keys()):
        P = f0.source.parent #Source of f
        Q = f0.target.parent #Target of f
    answer = SimpleDStructure(F2,P.algebra)        
    from_to_new = dict()
    to_to_new = dict()
    new_to_old = dict()
    #Add the generators
    for x in P.getGenerators():
        newx = DGenerator(answer, x.idem)
        from_to_new[x] = newx
        new_to_old[newx] = x
        answer.addGenerator(newx)
    for x in Q.getGenerators():
        newx = DGenerator(answer, x.idem)
        to_to_new[x] = newx
        new_to_old[newx] = x
        answer.addGenerator(newx)
    #Add the differential on P
    for x in P.getGenerators():
        dx = P.delta(x)
        for ay in list(dx.keys()):
            answer.addDelta(from_to_new[x],from_to_new[ay[1]],ay[0],dx[ay])
    for x in Q.getGenerators():
        dx = Q.delta(x)
        for ay in list(dx.keys()):
            answer.addDelta(to_to_new[x],to_to_new[ay[1]],ay[0],dx[ay])
    for f0 in list(f.keys()):
        answer.addDelta(from_to_new[f0.source],to_to_new[f0.target],f0.coeff,f[f0])
    if returnDicts:
        return (answer, from_to_new, to_to_new, new_to_old)
    return answer

def isQIDmor(f):
    "Check if a morphism of type D structures is a quasi-isomorphism."
    cone = mappingConeD(f)
    #Test if f is a chain map:
    if cone.testDelta():
        #Check if cone is acyclic
        cone.simplify()
        if len(cone) == 0:
            return True
    return False

def findQI(dMors):
    "Return an element of dMors which is a quasi-isomorphism."
    for f in dMors:
        if isQIDmor(f):
            return f
    assert False

def involutiveCx(P,Q, mult_one = True, sanityTests = False, verbose = False):
    """Returns the rank of involutive Floer homology of Y.s
    Input: P, Q: type D modules for handlebodies so that CF^(Y) = H_*(Mor(P,Q))
    """
    #How this works:
    #Compute: 
        #0. MP = azDA(pmc).tensorD(P), MQ = azDA(pmc).tensorD(Q)
        #1. A graded homotopy equivalence PhiPinv from P to MP
        #2. A graded homotopy equivalence PhiQ from MQ to Q
        #3. P.morToD(Q).simplify(find_homology_basis = True). Call result PQsimp
    #Let iota be the composition PQsimp -> (MP).morToD(MQ) -> P.morToD(Q)
    #given by f -> g = tensorDAid(P,Q,MP,MQ,f) -> PhiQ\circ g\circ PhiPinv
    #Compute iota on a basis of PQsimp. Result is a morphism of chain complexes.
    #Take mapping cone of Id+iota and take homology
    pmc = P.algebra.pmc
    M = azDA(pmc)
    if verbose:
        print("Lengths of (P,Q,M)"+repr((len(P.getGenerators()),len(Q.getGenerators()),len(M.getGenerators()))))
    MP = M.tensorD(P)
    if verbose:
        print("MP has length "+repr(len(MP.getGenerators())))
    MQ = M.tensorD(Q)
    if verbose:
        print("MQ has length "+repr(len(MQ.getGenerators())))
    PtoMPcx = P.morToD(MP) #This step tends to be very slow.
    if verbose:
        print("Computed PtoMPcx. Number of generators:"+repr(len(PtoMPcx.getGenerators())))
    PtoMQcx = P.morToD(MQ)
    if verbose:
        print("Computed PtoMQcx.  Number of generators:"+repr(len(PtoMQcx.getGenerators())))
    MQtoQcx = MQ.morToD(Q)
    if verbose:
        print("Computed MQtoQcx.  Number of generators:"+repr(len(MQtoQcx.getGenerators())))
    PtoMPcx.simplify(find_homology_basis = True) #This step also slow.
    if verbose:
        print("Simplified PtoMPcx. Number of generators: "+repr(len(PtoMPcx.getGenerators())))
    MQtoQcx.simplify(find_homology_basis = True)
    if verbose:
        print("Simplified MQtoQcx. Number of generators: "+repr(len(MQtoQcx.getGenerators())))

    PhiPinv = findQI([x.prev_meaning for x in PtoMPcx.getGenerators()])    
    if verbose:
        print("Found quasi-isomorphism PhiPinv")
    PhiQ = findQI([x.prev_meaning for x in MQtoQcx.getGenerators()])    
    if verbose:
        print("Found quasi-isomorphism PhiQ")
    if sanityTests:
        assert isQIDmor(PhiPinv)
        assert isQIDmor(PhiQ)
    PQ = P.morToD(Q)
    MPtoMQcx = MP.morToD(MQ) #Presumably this is the limiting step.
    if verbose:
        print("Computed PQ, MPtoMQcx. Number of generators of (PQ, MPtoMQcx):"+repr((len(PQ.getGenerators()),len(MPtoMQcx.getGenerators()))))
    PQ.simplify(find_homology_basis = True)
    if verbose:
        print("Simplified PQ. Number of generators: "+repr(len(PQ.getGenerators())))
    PQtarg = P.morToD(Q)
    iota = SimpleChainMorphism(PQ,PQtarg)    
    for gen_from in PQ.generators:
        gf_tens_id = tensorDAid(P,Q,M,MP,MQ,MPtoMQcx,gen_from.prev_meaning)
        if sanityTests:
            assert MPtoMQcx.diffElt(gf_tens_id) == E0
        comp1 = composeMor(gf_tens_id,PhiPinv, parent=PtoMQcx)
        if sanityTests:
            assert PtoMQcx.diffElt(comp1) == E0
        composedMap = composeMor(PhiQ,comp1, parent=PQtarg)
        if sanityTests:
            assert PQtarg.diffElt(composedMap) == E0
        for gen_to in list(composedMap.keys()):
            iota.addMorphism(gen_from, gen_to, composedMap[gen_to])
    if verbose:
        print("Computed iota")
    Id = SimpleChainMorphism(PQ,PQtarg)
    for gen_from in PQ.generators:
        for gen_to in list(gen_from.prev_meaning.keys()):
            new_gen_to = MorDtoDGenerator(PQtarg,gen_to.source,gen_to.coeff,gen_to.target)
            Id.addMorphism(gen_from, new_gen_to, gen_from.prev_meaning[gen_to])
    if sanityTests:
        #Sanity checks; should always pass
        assert isQIcx(Id)
        assert isQIcx(iota)
    idPlusIota = Id.sum(iota)
    if verbose:
        print("Computed Id + iota")
    mappingcone = idPlusIota.mappingConeCx()
    if verbose:
        print("Computed mapping cone.")
    mappingcone.simplify()
    return mappingcone

def invOfDCov(str_input, verbose=False):
    """Compute HF^ and HFI^ of the branched double cover of given knot.
    Input has the same form as readBridgePresentation (which can be found in data/input_12_FL.txt)"""
    bridgePres = readBridgePresentation(str_input)
    return invOfDCovSplit(str_input,len(bridgePres.braid_word),verbose)

def invOfDCovSplit(str_input, split_index, verbose=False):
    """Compute HF^ and HFI^ of the branched double cover of given knot, 
    splitting the diagram after split_index crossings to compute the involutive complex.
    Input has the same form as readBridgePresentation (which can be found in data/input_12_FL.txt)"""
    bridgePres = readBridgePresentation(str_input)
    braidGrp = Braid(bridgePres.num_strands)
    slides = braidGrp.getArcslides(bridgePres.braid_word)
    slidesInv = [x.inverse() for x in slides]
    slidesDA = [ArcslideDA(x) for x in slides]
    slidesInvDA = [ArcslideDA(x) for x in slidesInv]
    slidesDA.reverse()
    answer = list()
    cupD = BraidCap(bridgePres.start).openCap()
    capD = BraidCap(bridgePres.end).openCap()
    P = cupD
    Q = capD
    for slideDA in slidesInvDA[:split_index]:
        P = slideDA.tensorD(P)
        P.simplify()
    for slideDA in slidesDA[:len(slidesDA)-split_index]:
        Q = slideDA.tensorD(Q)
        Q.simplify()
    cx = P.morToD(Q)
    cx.simplify()
    if verbose:
        print("HF = "+repr(len(cx))+". Working on IHF next.")
    invcx = involutiveCx(P,Q)
    return invcx


def checkSplits(str_input):
    "Attempt to find best place to split up a bridge diagram for involutive computation."
    bridgePres = readBridgePresentation(str_input)
    braidGrp = Braid(bridgePres.num_strands)
    slides = braidGrp.getArcslides(bridgePres.braid_word)
    slidesInv = [x.inverse() for x in slides]
    slidesDA = [ArcslideDA(x) for x in slides]
    slidesInvDA = [ArcslideDA(x) for x in slidesInv]
    slidesDA.reverse()
    answer = list()
    cupD = BraidCap(bridgePres.start).openCap()
    capD = BraidCap(bridgePres.end).openCap()
    for i in range(len(bridgePres.braid_word)):
        P = cupD
        Q = capD
        for slideDA in slidesInvDA[:i]:
            P = slideDA.tensorD(P)
            P.simplify()
        for slideDA in slidesDA[:len(slidesDA)-i]:
            Q = slideDA.tensorD(Q)
            Q.simplify()
        cx = P.morToD(Q)
        cxlen = len(cx)
        cx.simplify()
        print((i,len(P.getGenerators()),len(Q.getGenerators()),cxlen,len(cx)))
        answer.append((i,len(P.getGenerators()),len(Q.getGenerators()),cxlen,len(cx)))
    return answer

