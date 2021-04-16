#!/usr/bin/env sage
from itertools import permutations

# TODO:
# Verification of design ness
# Affine Singer

def binom(n,r,q):
    '''
    Computes the Gaussian binomial coefficient (n r)_q. This happens to be the number of r dimensional subspaces of an n dimensional vector space over the finite field F_q.
    '''
    if n == r:
        return 1
    elif r < n:
        topfrac = 1
        for i in range(n-r+1,n+1):
            topfrac = topfrac * (1-q**i)
        bottomfrac=1
        for i in range(1,r+1):
            bottomfrac = bottomfrac * (1-q**i)
        bino = float(topfrac/bottomfrac)
        assert(bino.is_integer()), "binom returns an int"
        return int(bino)

class designblock():
    '''
    A block of a design
    '''
    # This was done so that changing blocks updates things like the resolution classes
    def __init__(self, elements):
        self.elements = elements

    def setminus(self, block):
        return designblock([i for i in self.elements if i not in block.elements])

    def intersection(self, block):
        return designblock([i for i in self.elements if i in block.elements])

class BIBD():
    '''
    Balanced incomplete block design. A balanced incomplete block design is formed from a set of points V
    and a collection of k-subsets of V, called B, with the property that every two points in V occur in the
    same number of elements of B. For instance, if V={0,1,2,3,4,5,6} and B is given by
    {{1,2,4},
    {2,3,5},
    {3,4,6},
    {4,5,0},
    {5,6,1},
    {6,0,2},
    {0,1,3}}
    has any two elements of V (for instance, 2 and 5) occurring in exactly one block ({2,3,5}).
    '''
    def __init__(self, parameters):
        if parameters != []:
            self.parameters = parameters
            self.v, self.k, self.lambduh = (int(parameters[0]), int(parameters[1]), int(parameters[2]))
            # v is the size of V, k is how long each block is, and lambduh (lambda) is how many blocks each pair occurs in
            self.calculate_extra_parameters()
        

    def calculate_extra_parameters(self):
        self.parameters = (self.v, self.k, self.lambduh)
        self.b = self.lambduh * self.v * (self.v-1)/(self.k*(self.k-1))
        self.r = self.lambduh * (self.v-1)/(self.k-1)
        # b is the number of blocks, and r is how many times each point occurs (it is the same for each point)
        assert(self.b.is_integer() and self.r.is_integer()), "inadmissible design"
        if self.b == self.v:
            self.symmetric = True
        else:
            self.symmetric = False

    def derived_params(self):
        '''
        Calculates the parameters for the derived design. A derived design is the intersection of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.symmetric == True and self.lambduh > 1:
            return [self.k, self.lambduh, self.lambduh - 1]
    
    def derived(self, blockstar=-1):
        '''
        Returns the residual design with respect to blockstar. Defaults to final block.
        '''
        if blockstar == -1:
            blockstar = self.blocks[-1]
        derived_design = BIBD(self.derived_params())
        iteratelist = self.blocks
        iteratelist.remove(blockstar)
        derived_design.blocks = [block.intersection(blockstar) for block in iteratelist]
        return derived_design

    def residual_params(self):
        '''
        Calculates the parameters for the residual design. The residual design is the difference of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.symmetric == True:
            return [self.v-self.k, self.k - self.lambduh, self.lambduh]
    
    def residual(self, blockstar=-1):
        '''
        Returns the residual design with respect to blockstar. Defaults to final block.
        '''
        if blockstar == -1:
            blockstar = self.blocks[-1]
        residual_design = BIBD(self.residual_params())
        iteratelist = self.blocks
        iteratelist.remove(blockstar)
        residual_design.blocks = [block.setminus(blockstar) for block in iteratelist]
        return residual_design

    def complement_params(self):
        '''
        Calculates the parameters for the complement design. The complement design is the complement of 
        all of the blocks in the total point set.
        '''
        return [self.v, self.v-self.k, self.b-2*self.r+self.lambduh]

    def complement(self, blockstar=-1):
        '''
        Returns the residual design with respect to blockstar. Defaults to final block.
        '''
        if blockstar == -1:
            blockstar = self.blocks[-1]
        complement_design = BIBD(self.complement_params())
        iteratelist = self.blocks
        iteratelist.remove(blockstar)

        points = []
        for block in self.blocks:
            for point in block.elements:
                if point not in points:
                    points += [point]
        pointset = designblock(points)
        complement_design.blocks = [pointset.setminus(block) for block in iteratelist]
        return complement_design

    def biject_obvious(self):
        '''
        Finds a bijection from the current set of points to a subset of the integers
        '''
        points = []
        for block in self.blocks:
            for point in block.elements:
                if point not in points:
                    points += [point]
        for block in self.blocks:
            block.elements = [points.index(point) for point in block.elements]

    def list_blocks(self):
        '''
            Returns the blocks of a design, in list form
        '''
        return [block.elements for block in self.blocks]

    def list_resolution_classes(self):
        '''
            Returns the resolution classes of a design, in list form
        '''
        return [[block.elements for block in resolclass] for resolclass in resolutionclasses]

class AG(BIBD):
    '''
    An affine geometry of order n over F_q. An affine geometry is the set of all d-flats of an n dimensional vector space over F_q. And a d-flat is a coset of a d-dimensional subspace. An example of an infinite affine geometry (which this program is not concerned with) is the set of all lines in a 3 dimensional space. Or the set of all planes. A line or plane through the origin is a 1 or 2 dimensional subspace, and a general line or plane is a translated version of that, so a coset. The points of the affine geometry are just the 0-flats.
    '''
    def __init__(self,  n, q, d=-1):
        if d == -1:
            d = n-1
        self.n = n
        self.q = q
        self.d = d
        self.v, self.k, self.lambduh = (q**n,   q**d,  binom(n-1, d-1, q))
        self.parameters = [int(self.v), int(self.k), int(self.lambduh)]
        BIBD.__init__(self, self.parameters)

    def getcosets(self, subspace, vectorspace):
        '''
        Gets the cosets of a given subspace in a vectorspace.
        '''
        availpoints = [i for i in vectorspace]
        cosets = []
        while availpoints != []:
            point = availpoints[0]
            coset = [point + elem for elem in subspace]
            cosets += [coset]
            [availpoints.remove(element) for element in coset]
        return cosets

    def bijection(self, point, q):
        '''
        Converts a tuple of length n over F_q into an integer, with the obvious bijection
        '''
        output = 1
        if type(point[0]) == sage.rings.finite_rings.integer_mod.IntegerMod_int:
            for i in range(len(point)):
                output += q**i * int(point[i])
            return output
        elif type(point[0]) == sage.rings.finite_rings.element_givaro.FiniteField_givaroElement:
            for i in range(len(point)):
                output += q**i * point[i].integer_representation()
            return output

    def biject(self):
        '''
        Bijects the affine geometry to the integers, mapping any tuple (a_{n-1}, a_{n-2}, ... a_0) to sum q^i a_i.
        If a_i are polynomials, converts them to integers.
        '''
        for block in self.blocks:
            block.elements = [self.bijection(point,self.q) for point in block.elements]

    def generate(self):
        '''
        Generates the affine geometry explicitly
        '''

        field = GF(self.q)
        V = VectorSpace(field, self.n)

        self.resolutionclasses = [[designblock(coset) for coset in self.getcosets([point for point in subspace], V)] for subspace in V.subspaces(self.d)]

        blocks = []
        for resolutionclass in self.resolutionclasses:
            blocks += resolutionclass
        self.blocks = blocks


class PG(BIBD):
    '''
    A projective geometry of order n over F_q. Given an n+1 dimensional vector space V over F_q, the projection of a k+1 dimensional subspace U is defined
    to be the set of all one dimensional subspaces of U. This has dimension k. Alternatively, we could define an equivalence relation ~ that is "being a scalar 
    multiple of each other", and the projection of U would be U/~. The projective geometry of order n over F_q has its k dimensional subspaces given by the 
    projections of k+1 dimensional subspaces in an n+1 dimensional vector space over F_q. If you don't know what finite projective geometries are, I honestly 
    don't know if this would help. 
    '''
    def __init__(self,  n, q, d=-1):
        if d == -1:
            d = n-1
        self.n, self.q, self.d = n, q, d
        self.v, self.k, self.lambduh = (binom(n+1,  1, q), binom(n, 1, q), binom(n-1, 1, q))
        self.parameters = (int(self.v), int(self.k), int(self.lambduh))
        BIBD.__init__(self, self.parameters)

    def projection(self, U, space):
        '''
        Finds the projection of U in the space
        '''
        return [space.subspace([point]) for point in U]

    def generate(self):
        '''
        Generates the projective geometry explicitly
        Note: blocks = [...] creates q(^k?) too many of each element. Optimisation needed.
        '''
        field = GF(self.q)
        V = VectorSpace(field, self.n+1)
        
        subspaces = [subspace for subspace in V.subspaces(self.d+1)]
        projections = [self.projection([i for i in subspace], V) for subspace in V.subspaces(self.d+1)]
        blocks = []
        for projection in projections:
            block = []
            for space in projection:
                if space not in block:
                    block += [space]
            blocks += [block]
        self.blocks = [designblock([[elem for elem in space] for space in block if space.dimension() != 0]) for block in blocks]

    def generate_hyperplanes(self):
        '''
        Generates the design isomorphic to the hyperplanes of the projective geometry
        '''
        # Uses Singer's method. Since there is an identification from an n+1 dimensional vector space over F_q to F_{q^(n+1)}, and
        # multiplication by x in F_{q^(n+1)} maps subspaces to subspaces, it is an automorphism of the design. Since we can identify
        # each element of the finite field Q(x) with a number such that x^i = Q(x), we can create an additive difference set to generate
        # our geometry with.

        field = GF(self.q**(self.n+1))
        powers = field.primitive_element().powers(self.q**(self.n+1))
        p = field.base().order()
        alpha = log(self.q)/log(p)
        diffset = [i for i in range(len(powers)) if powers[i].polynomial().degree() <= self.n*alpha -1 and i < (self.q**(self.n+1)-1)/(self.q-1)]
        diffdesign = difference_method([diffset], (self.q**(self.n+1)-1)/(self.q-1))
        diffdesign.generate()
        self.blocks = diffdesign.blocks

class difference_method(BIBD):
    def __init__(self, difference_sets, n, repetition="False"):
        self.difference_sets = difference_sets
        self.n = n
        self.repetition = repetition
        self.verify()
        
    def verify(self):
        differences = [0 for i in range(self.n)]
        for diffset in self.difference_sets:
            for x,y in permutations(diffset, 2):
                try:
                    differences[(x-y)%self.n] += 1
                except KeyError:
                    differences[(x-y)%self.n] = 1
        differences[0]=differences[1]

        assert(len(set(differences)) == 1), "Not a difference system: Unequal differences"
        assert(len(set([len(i) for i in self.difference_sets])) == 1), "Not a difference system: Lengths not fixed"

        self.lambduh = differences[0]
        self.k = len(self.difference_sets[0])
        self.parameters = (self.n, self.k, self.lambduh)
        BIBD.__init__(self, self.parameters)
            

    def generate(self):
        sets = []
        for difference_set in self.difference_sets:
            sets += [[(x+i)%self.n for x in difference_set] for i in range(self.n)]
        if self.repetition == False:
            self.blocks = [designblock(theset) for theset in list(set(sets))]
        else:
            self.blocks = [designblock(theset) for theset in sets]

def quad_residue_diff_set(q):
    if q%4 == 3:
        pass
    else:
        print("Error: q != 3 mod 4")
    residues = quadratic_residues(q)
    residues.remove(0)
    return residues

def quad_residue_design(q):
    return difference_method([quad_residue_diff_set(q)],q)
