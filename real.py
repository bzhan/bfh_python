from pmc import PMC, Idempotent, connectSumPMC, StrandAlgebra, Strands, StrandDiagram, splitPMC
from dstructure import SimpleDStructure, SimpleDGenerator
from dastructure import SimpleDAStructure, SimpleDAGenerator
from algebra import F2
from ddstructure import identityDD

def realAZ(pmc: PMC) -> SimpleDStructure:
    "Returns the small model for the real Auroux-Zarev piece of the double pmc#(-pmc)."
    # Construct the pointed matched circle for the boundary
    N: int = 2*pmc.n
    pmc_alg: StrandAlgebra = pmc.getAlgebra(mult_one = True)
    pmc2: PMC = connectSumPMC(pmc,pmc.opp())
    algebra: StrandAlgebra = pmc2.getAlgebra()
    dstr = SimpleDStructure(F2, algebra)

    alg_to_gens = dict()

    # Construct the generators
    for a in pmc.getStrandDiagrams(pmc_alg):
        # Get left idempotent. This is complement of left idem of a union r of that.
        aidem_left: Idempotent = a.left_idem.comp() #This is a list of pairids
        aidem_right: Idempotent = a.right_idem.comp() #This is a list of pairids
        pairids = list()
        for i in aidem_left:
            pmcpair: list[tuple[int]] = pmc.pairs[i]
            pairids.append(pmc2.pairid[pmcpair[0]])
        for i in aidem_right:
            pmcpair = pmc.pairs[i]
            pairids.append(pmc2.pairid[N - pmcpair[0]-1])
        idem = Idempotent(pmc2,pairids)
        gen = SimpleDGenerator(dstr, idem, repr(a))
        alg_to_gens[a] = gen
        dstr.addGenerator(gen)

    # Type (ii) differentials
    a: StrandDiagram
    for a in pmc.getStrandDiagrams(pmc_alg):
        b = pmc_alg.diff(a)
        c: StrandDiagram
        for c in b.keys():
            dstr.addDelta(alg_to_gens[a],alg_to_gens[c],alg_to_gens[a].idem.toAlgElt(algebra),1)

    # Build a list of chords
    chords = [Strands(pmc, [(i, j)]) for i in range(pmc.n) for j in range(i+1, pmc.n)]

    # Type (vi) differentials
    for a in pmc.getStrandDiagrams(pmc_alg):
        for xi in chords:
            mirror_strands = Strands(pmc2, [(N-j-1,N-i-1) for (i,j) in xi])
            if mirror_strands.leftCompatible(alg_to_gens[a].idem) and xi.leftCompatible(a.right_idem):
                xialg = StrandDiagram(pmc_alg,a.right_idem,xi)
                b: StrandDiagram
                for b in a*xialg:
                    xi_alg_elt = StrandDiagram(algebra, alg_to_gens[a].idem, mirror_strands)
                    if alg_to_gens[b].idem == xi_alg_elt.right_idem:
                        dstr.addDelta(alg_to_gens[a],alg_to_gens[b],xi_alg_elt,1)

    # Type (vii) differentials
    for a in pmc.getStrandDiagrams(pmc_alg):
        for xi in chords:
            these_strands = Strands(pmc2, [(i,j) for (i,j) in xi])
            mirror_strands = Strands(pmc2, [(N-j-1,N-i-1) for (i,j) in xi])
            if these_strands.leftCompatible(alg_to_gens[a].idem) and xi.rightCompatible(a.left_idem):
                xialg = StrandDiagram(pmc_alg,left_idem=None,strands=xi,right_idem=a.left_idem)
                b: StrandDiagram
                for b in xialg*a:
                    xi_alg_elt = StrandDiagram(algebra, alg_to_gens[a].idem, these_strands)
                    if alg_to_gens[b].idem == xi_alg_elt.right_idem:
                        dstr.addDelta(alg_to_gens[a],alg_to_gens[b],xi_alg_elt,1)
    return dstr

def realAZda(pmc: PMC) -> SimpleDAStructure:
    "Returns the type A structure for the small model of the real Auroux-Zarev piece of the double pmc#(-pmc). Returns a da structure where left algebra is splitPMC(0), because type A modules haven't been implemented on their own."
    pmc_alg: StrandAlgebra = pmc.getAlgebra(mult_one = True)
    pmc2: PMC = connectSumPMC(pmc,pmc.opp())
    algebra: StrandAlgebra = pmc2.getAlgebra()
    trivial_alg: StrandAlgebra = splitPMC(0).getAlgebra()
    one = list(trivial_alg.getGenerators())[0]
    one_idem = list(trivial_alg.getIdempotents())[0]
    astr = SimpleDAStructure(F2, trivial_alg, algebra)
    alg_to_gens = dict()
    #Add the generators
    for a in pmc_alg.getGenerators():
        a: StrandDiagram
        idem = connectSumIdem(a.right_idem,a.left_idem.opp(),pmc2)
        gen = SimpleDAGenerator(astr, one_idem, idem, repr(a))
        astr.addGenerator(gen)
        alg_to_gens[a] = gen
    #Add the operations
    #Differentials from smoothing a crossing in a:
    for a in pmc_alg.getGenerators():
        a: StrandDiagram
        for b in a.diff():
            astr.addDelta(alg_to_gens[a],alg_to_gens[b],one,[],1)
    #m_2 operations, coming from splitting up an algebra element and multiplying
    for a in pmc_alg.getGenerators():
        for b in pmc_alg.getGenerators():
            for c in pmc_alg.getGenerators():
                if(a*b*c and not (a.isIdempotent() and c.isIdempotent()) ):
                    for d in a*b*c:
                        #Turn (a,c) into an element of the bigger algebra
                        #Then generator corresponding to b * generator corresponding to (a,c) is generator corresponding to d
                        start_gen = alg_to_gens[b]
                        end_gen = alg_to_gens[d]
                        e = connectSumStrandDiagram(c,a.opp(),algebra)
                        assert (e.left_idem == start_gen.idem2 and e.right_idem == end_gen.idem2)
                        astr.addDelta(start_gen,end_gen,one,[e],1)
    return astr


def connectSumIdem(idem1:Idempotent,idem2:Idempotent,parent:PMC) -> Idempotent:
    "Give idempotent in connected sum PMC parent from taking connected sum of these idempotents."
    idem = list(idem1) + [p+idem1.pmc.num_pair for p in idem2]
    return Idempotent(parent, idem)

def connectSumStrandDiagram(a:StrandDiagram, b:StrandDiagram, parent:StrandAlgebra) -> StrandDiagram:
    "Give strand diagram in StrandAlgebra parent (algebra of connected sum of PMCs for a and b) from taking connected sum of a and b"
    n = a.pmc.n
    strands = Strands(parent.pmc, tuple(a.strands)+tuple((i+n,j+n) for (i,j) in b.strands))
    return StrandDiagram(parent, connectSumIdem(a.left_idem,b.left_idem,parent.pmc),strands)

def selfGluingHB(pmc: PMC) -> SimpleDStructure:
    "Currently hardwires multiplicity 1"
    new_pmc: PMC = connectSumPMC(pmc,pmc.opp()) 
    new_alg: StrandAlgebra = new_pmc.getAlgebra(mult_one=True)
    answer = SimpleDStructure(F2,new_alg)
    # Will import most of the data from dd identity Need them in all sizes of idempotent
    dd_gens_to_new_gens = dict()
    idems_to_new_gens = dict() #The idempotent also determines the generator
    for idem_size in range(0,2*pmc.genus+1):
        dd_id = identityDD(pmc,idem_size)
        # Import the generators
        for x in dd_id.generators:
            newx = SimpleDGenerator(answer, connectSumIdem(x.idem1,x.idem2,new_pmc),(idem_size,x.name))
            dd_gens_to_new_gens[x] = newx
            idems_to_new_gens[newx.idem] = newx
            answer.addGenerator(newx)
        # Import the differentials
        for x in dd_id.generators:
            for (a,b,y) in dd_id.delta(x):
                c = connectSumStrandDiagram(a,b,new_alg)
                answer.addDelta(dd_gens_to_new_gens[x],dd_gens_to_new_gens[y],c,1) #Seems okay to here.
    # Add the symmetric chords
    chords = [Strands(new_pmc,[(i,new_pmc.n-i-1),]) for i in range(pmc.n)]
    for x in answer.generators:
        for xi in chords:
            if xi.leftCompatible(x.idem):
                alg_coeff = StrandDiagram(new_alg,x.idem,xi)
                right_idem = alg_coeff.right_idem
                y = idems_to_new_gens[right_idem]
                answer.addDelta(x,y,alg_coeff,1)
    return answer