##!/usr/bin/env sage
#from itertools import permutations
#
# TODO:
# Verification of design ness
# Affine Singer
# Documentation of sage classes and input/output/example
# Documentation of functions
# More example usage in README
# Extend affine to projective, rename projective to affine
# Speed up PG.generate() and other timing based stuff
# Fix PG.generate_hyperplanes for prime power
# Generate designs that exist

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

class designpoint():
    '''
    A point of the design
    '''
    # This was done so that bijections automatically update things
    def __init__(self, value):
        self.value = value

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

class existinfo():
    '''
    Wrapper to hold existence information
    '''
    def __init__(self, exist=0.5, method = "", parameters=[], message=""):
        self.exist = exist
        self.method = method
        self.parameters = parameters
        self.message = message

class BIBD():
    '''
    Balanced incomplete block design. A balanced incomplete block design is formed from a set of points V
    and a collection of k-subsets of V, called B, with the property that every two points in V occur in the
    same number of elements of B. For instance, if V={0,1,2,3,4,5,6} and B is given by
    {{1,2,4},
    {)2,3,5},
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
            self.is_complement = False
            self.is_residual = False
            self.is_residual_inverse = False
        

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

    def pointinV(self, point):
        '''
        Returns the point object in the set of the designs points with the given value.
        '''
        for v in self.V:
            if point == v.value:
                return v

    ### Creating new designs ###
    def derived_params(self):
        '''
        Calculates the parameters for the derived design. A derived design is the intersection of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.symmetric == True and self.lambduh > 1:
            return [self.k, self.lambduh, self.lambduh - 1]

    def derived_inverse_params(self):
        x = self.complement()
        x = x.residual()
        x = x.complement()
        return x.parameters

    def derived(self, blockstar=-1):
        return derived_design(self, blockstar)

    def derived_inverse(self):
        return BIBD([derived_inverse_params(self)])

    def residual_params(self):
        '''
        Calculates the parameters for the residual design. The residual design is the difference of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.symmetric == True:
            return [self.v-self.k, self.k - self.lambduh, self.lambduh]

    def residual_inverse_params(self):
        if self.r == self.k + self.lambduh:
            return [self.v+self.k+self.lambduh, self.k + self.lambduh, self.lambduh]

    def residual_inverse(self):
        return BIBD(self.residual_inverse_params())

    def residual(self, blockstar=-1):
        return residual_design(self, blockstar)
    
    def complement_params(self):
        '''
        Calculates the parameters for the complement design. The complement design is the complement of 
        all of the blocks in the total point set.
        '''
        return [self.v, self.v-self.k, self.b-2*self.r+self.lambduh]

    def complement(self):
        return complement_design(self)

    def multiple(self, n):
        return multiple_design(self, n)

    ### Bijection ###
    def biject_obvious(self):
        '''
        Finds a bijection from the current set of points to a subset of the integers
        '''
        for i in range(len(self.V)):
            self.V[i].value = i
    
    def decouple_points(self):
        self.V = [designpoint(point.value) for point in self.parent.V]
        for block in self.blocks:
            for i in range(len(block.elements)):
                block.elements[i] = self.pointinV(block.elements[i].value)

    ### Listing information ###
    def list_points(self):
        '''
            Returns the points of a design, in list form
        '''
        return [point.value for point in self.V]

    def list_blocks(self):
        '''
            Returns the blocks of a design, in list form
        '''
        return [[point.value for point in block.elements] for block in self.blocks]

    def list_resolution_classes(self):
        '''
            Returns the resolution classes of a design, in list form
        '''
        return [[[point.value for point in block.elements] for block in resolclass] for resolclass in self.resolutionclasses]


    ### Basic existence methods ###
    def exists_Hadamard_matrix(self, n):
        if n == 0:
            exist = existinfo(exist=False, message="Trivially")
            return exist

        if n == 2:
            exist = existinfo(exist=True, message="Hadamard matrix of order 2 exists")
            return exist

        if n%4 != 0:
            exist = existinfo(exist=False, message = "Hadamard matrix order not divisible by 4")
            return exist

        if is_prime_power(n-1) and (n-1)%4 == 3:
            exist = existinfo(exist=True, message = str(n-1) + " is prime power 3 mod 4")
            return exist

        if n%2 == 0:
            exist = self.exists_Hadamard_matrix(n/2)
            if exist.exist == True:
                exist2 = existinfo(exist=True, message = str("times 2, ") + exist.message)
                return exist2
        
        if n%4 == 0:
            for divisor1 in divisors(n):
                divisor2 = n/(divisor1)
                exist1 = self.exists_Hadamard_matrix(divisor1)
                exist2 = self.exists_Hadamard_matrix(divisor2)
                if exist1.exist == True and exist2.exist == True:
                    exist = existinfo(exist=True, message = str("Formed from: ") + exist1[1] + ", and " + exist2[1])
                    return exist

        exist = existinfo(message = "Tried everything")
        return exist

    def permits_Hadamard_design(self):
        n = self.k - self.lambduh
        exist = self.exists_Hadamard_matrix(4*n)
        if exist.exist == True:
            if self.v == 4*n-1 and self.k == 2*n-1 and self.lambduh == n-1:
                exist2 = existinfo(exist=True, message=exist.message, parameters = [n])
                return exist2
        return existinfo(exist=False)

    def permits_AG(self):
        if is_prime_power(self.v):
            p1, alpha1 = is_prime_power(self.v, get_data=True)
        else:
            return existinfo(exist=False)

        if is_prime_power(self.k):
            p2, alpha2 = is_prime_power(self.k, get_data=True)
        else:
            return existinfo(exist=False)

        if p1 != p2:
            return existinfo(exist=False)
        maximalq = gcd(alpha1, alpha2)
        for divisor in maximalq.divisors():
            if binom(alpha1/divisor - 1, alpha2/divisor-1, p1**divisor) == self.lambduh:
                return existinfo(exist=True, parameters = [alpha1/divisor, p1**divisor, alpha2/divisor], method=AG)
        return existinfo(exist=False)

    def permits_PG(self):
        '''
        Optimise
        '''
        possiblesols = []
        primepower = 2
        n = 1
        cumulativeresults = []
        while binom(2,1,primepower) <= self.v:
            guessselfn = binom(n,1,primepower)
            if guessselfn == self.v:
                possiblesols += [[n-1,primepower, cumulativeresults]]
                primepower = next_prime_power(primepower)
                n = 1 
                cumulativeresults = []
                
            if guessselfn < self.v:
                n += 1
                cumulativeresults += [guessselfn]
            else:
                primepower = next_prime_power(primepower)
                n = 1 
                cumulativeresults = []
        
        for n, primepower, cumulativeresults in possiblesols:
            if self.k in cumulativeresults:
                if binom(n-1,cumulativeresults.index(self.k)-1,primepower) == self.lambduh:
                    return existinfo(exist=True, parameters = [n,primepower,cumulativeresults.index(self.k)], method=PG)
        return existinfo(exist=False)
         
    def permits_quad_residue(self):
        if is_prime_power(self.v):
            if self.k == (self.v-1)/2 and self.lambduh == (self.v-3)/4:
                return existinfo(exist=True, parameters = [self.v], method = quad_residue_design)

        return existinfo(exist = False)
    
    ### Complicated existence methods
    def existence(self):
        '''
        Checks if a design can be constructed as a constructible complement or derived or multiple design
        '''

        basicexistence = self.existence_basic()
        if basicexistence.exist == True:
            return existinfo(exist=True, message = basicexistence.message)

        if basicexistence.exist == False:
            return existinfo(False, message = basicexistence.message)

        for possiblelambduh in divisors(self.lambduh):
            if possiblelambduh != self.lambduh:
                try:
                    x = BIBD([self.v,self.k,possiblelambduh])
                    exist = x.existence()
                    if exist.exist == True:
                        return existinfo(exist=True, message = "Multiply by " + str(self.lambduh/possiblelambduh) + ": " + exist.message)
                except AssertionError:
                    pass

        if self.is_complement == False and self.v - self.k > 1:
            y = self.complement()
            y.is_complement = True
            exist = y.existence()
            if exist.exist == True:
                return existinfo(exist=True, message = "Complement of: " + exist.message)

        if self.r == self.k + self.lambduh and self.is_residual == False:
            y = self.residual_inverse()
            y.is_residual_inverse = True
            exist = y.existence()
            if exist.exist == True:
                return existinfo(exist=True, message = "Residual of: " + exist.message)
            if self.lambduh < 3 and exist.exist == False:
                return existinfo(exist=False, message = "Quasiresidual designs of lambda < 3 are residual designs, but no residual exists: " + exist.message)

        if self.symmetric and self.is_residual_inverse == False and self.k - self.lambduh > 1:
            x = self.residual()
            x.is_residual = True
            if x.lambduh < 3:
                exist = x.existence()
                if exist.exist == True:
                    return existinfo(exist=True, message = "Quasiresidual design: " + str(x.parameters) + exist.message)

        return existinfo(exist=0.5, message = "Maybe")

    def existence_basic(self):
        '''
        Checks if a design can be constructed from other objects - vector fields, Hadamard matrices, etc
        '''
        a = self.k - self.lambduh
        b = (-1)**((self.v-1)/2)* self.lambduh

        hadamard = self.permits_Hadamard_design()
        if self.b < self.v:
            return existinfo(exist = False, message = "Violates Fischer's Inequality")

        elif self.k == self.v - 1 and self.lambduh == self.v - 2:
            return existinfo(exist=True, message = "Trivial case")

        elif self.permits_quad_residue().exist:
            return existinfo(exist=True, message = "Quadratic Residue " + str(self.permits_quad_residue().parameters))

        elif self.permits_AG().exist:
            return existinfo(exist=True, message = "Affine geometry " + str(self.permits_AG().parameters))

        elif self.permits_PG().exist:
            return existinfo(exist=True, message = "Projective geometry " + str(self.permits_PG().parameters))

        elif hadamard.exist:
            return existinfo(exist=True, message = "Hadamard design of order " + str(hadamard[0]) + ": Formed by Hadamard matrix of order " + str(4*hadamard[0]) + ": " + hadamard[1])

        if self.symmetric == True:
            if self.v%2 == 0:
                # The Bruck-Ryser-Chowla theorem states that, for v even, a (v,k,lambda) design exists, then
                # (k-lambda) is a perfect square
                if sqrt(a).is_integer():
                    return existinfo(exist=0.5, message = str(self.k-self.lambduh) + " is a perfect square")
                else:
                    return existinfo(exist=False, message = "by BRC, " + str(self.k-self.lambduh) + " not perfect square")

            elif self.v%2 == 1:
                # The Bruck-Ryser-Chowla theorem states that if, for v odd, a (v,k,lambda) design exists, then
                # z^2 = (k-lambda) x^2 + (-1)^((v-1)/2) lambda y^2 has a nontrivial solution


                # Let m be a factor appearing once in ab. Let c be b if m divides a, and a if b divides m. 
                # If {x^2 | x in Z_m} cap {c*x^2 | x in Z_m} = {0}, then the equation z^2 = ax^2 + by^2
                # has no solutions. Simply take both sides mod m, deduce two squares must be 0, and sub in.
                if a*b != 0:
                    possiblem = [fact[0] for fact in factor(a*b) if fact[1]==1]
                    for m in possiblem:
                        if m.divides(a):
                            c = int(b)
                        else:
                            c = int(a)

                        field = Integers(m)
                        possiblez = set([i**2 for i in field]).intersection(set([c*i**2 for i in field]))
                        if possiblez == {0}:
                            return existinfo(exist=False, message = "common divisibility argument mod " + str(m))
                else:
                    return existinfo(exist=0.5, message = "solution to BRC given by " + str((0,1,0)) + " or " + str((0,0,1)))

                # Easy solutions
                if sqrt(a).is_integer():
                    return existinfo(exist=0.5, message = "solution to BRC given by " + str((sqrt(a), 1, 0)))
                elif sqrt(self.lambduh).is_integer() and b>0:
                    return existinfo(exist=0.5, message = "solution to BRC given by " + str((sqrt(self.lambduh), 0, 1)))
                elif sqrt(a+b).is_integer():
                    return existinfo(exist=0.5, message = "solution to BRC given by " + str((sqrt(a+b), 1, 1)))
                elif a == 1 and b == 1:
                    return existinfo(exist=0.5, message = "solution to BRC given by " + str((5,3,4)))

                # If a=b=0 mod m and a/m = b/m != 0 mod m, we immediately know that z^2 = 0 mod m
                # and then a/m x^2 + b/m y^2 = 0 mod m, which means a/m (x^2+y^2) = 0 mod m. For m=3,
                # then x^2+y^2 in {0,1,2}, and since {x in Z_3 | exists c in Z_3\{0} : cx = 0} = {0},
                # we know that x^2+y^2 in {0}, so we have a common divisibility argument.
                elif a%3 == 0 and b%3 == 0 and (a/3)%3 == (b/3)%3 and ((a/3)%3 == 2 or (a/3)%3 == 1):
                    return existinfo(exist=False, message = "common divisibility argument mod 3")
                elif (a%4 == 0 and (a/4)%4 != 0 and (b%4 == 2 or b%4 == 3)) or (b%4 == 0 and (b/4)%4 != 0 and (a%4 == 2 or a%4 == 3)):
                    return existinfo(exist=False, message = "common divisibility argument mod 4")

                # Giving up
                else:
                    return existinfo(exist=0.5, message = "try BRC?" + str(self.parameters) + "z^2 = " + str(a) + "x^2 + " + str(b) + "y^2")
        else:
            return existinfo(exist=0.5, message = "maybe?")

### Classes for creating new designs from old ones ###
class derived_design(BIBD):
    def __init__(self, parent, blockstar):
        self.parent = parent
        self.blockstar = blockstar
        self.parameters = self.derived_params(parent)
        BIBD.__init__(self, self.parameters)

    def derived_params(self, parent):
        '''
        Calculates the parameters for the derived design. A derived design is the intersection of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.parent.symmetric == True and self.parent.lambduh > 1:
            return [self.parent.k, self.parent.lambduh, self.parent.lambduh - 1]
        else: 
            return []
    
    def generate(self, blockstar=-1):
        '''
        Returns the derived design with respect to blockstar. Defaults to final block.
        '''
        if blockstar == -1:
            blockstar = self.parent.blocks[-1]
        iteratelist = [i for i in self.parent.blocks]
        iteratelist.remove(blockstar)
        self.V = [point for point in self.parent.V if point in blockstar.elements]
        self.blocks = [block.intersection(blockstar) for block in iteratelist]

class residual_design(BIBD):
    def __init__(self, parent, blockstar):
        self.parent = parent
        self.blockstar = blockstar
        self.parameters = self.residual_params()
        BIBD.__init__(self, self.parameters)
        self.is_residual = True

    def residual_params(self):
        '''
        Calculates the parameters for the residual design. The residual design is the difference of all
        blocks with respect to some fixed block of the design. Only a design if symmetric.
        '''
        if self.parent.symmetric == True:
            return [self.parent.v-self.parent.k, self.parent.k - self.parent.lambduh, self.parent.lambduh]
        else:
            return []

    def generate(self, blockstar=-1):
        '''
        Returns the residual design with respect to blockstar. Defaults to final block.
        '''
        if blockstar == -1:
            blockstar = self.parent.blocks[-1]
        params = self.residual_params()
        assert(params != None), "Not symmetric!"
        iteratelist = [i for i in self.parent.blocks]
        iteratelist.remove(blockstar)
        self.V = [point for point in self.parent.V if point not in blockstar.elements]
        self.blocks = [block.setminus(blockstar) for block in iteratelist]

class complement_design(BIBD):
    def __init__(self, parent):
        self.parent = parent
        self.parameters = self.complement_params()
        BIBD.__init__(self, self.parameters)
        self.is_complement = True

    def complement_params(self):
        '''
        Calculates the parameters for the complement design. The complement design is the complement of 
        all of the blocks in the total point set.
        '''
        return [self.parent.v, self.parent.v-self.parent.k, self.parent.b-2*self.parent.r+self.parent.lambduh]

    def generate(self, blockstar=-1):
        '''
        Returns the complement design with respect to blockstar. 
        '''
        if blockstar == -1:
            blockstar = self.parent.blocks[-1]
        iteratelist = [i for i in self.parent.blocks]

        points = self.parent.V
        pointset = designblock(points)
        self.V = [i for i in points]
        self.blocks = [pointset.setminus(block) for block in iteratelist]

class multiple_design(BIBD):
    def __init__(self, parent, n):
        self.parent = parent
        self.parameters = [parent.v, parent.k, parent.lambduh * n]
        self.n = n
        BIBD.__init__(self, self.parameters)

    def generate(self):
        self.V = self.parent.V
        self.blocks = [i for i in self.parent.blocks]
        self.blocks = self.blocks * self.n

### Classes for creating new designs ###
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
        availpoints = [i for i in self.V]
        cosets = []
        while availpoints != []:
            point = availpoints[0].value
            coset = [self.pointinV(point + elem) for elem in subspace]
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
        for point in self.V:
            point.value = self.bijection(point.value, self.q)

    def generate(self):
        '''
        Generates the affine geometry explicitly
        '''

        field = GF(self.q)
        V = VectorSpace(field, self.n)
        self.V = [designpoint(v) for v in V]

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

    def projection_points(self, space):
        self.V = [designpoint(i) for i in space.subspaces(1)]

    def projection(self, U):
        '''
        Finds the projection of U in the space self.V
        '''
        projection = []
        # Note that what follows isn't a duplication of pointinV: the point in U is not a point
        # of the design, but a point of the vector space that is contained in a point of the design
        for point in U:
            for v in self.V:
                if point in v.value:
                    if v not in projection:
                        projection += [v]
                    break
        return designblock(projection)

    def generate(self):
        '''
        Generates the projective geometry explicitly
        Note: projections = [...] creates q(^k?) too many of each element. Optimisation needed.
        '''
        field = GF(self.q)
        V = VectorSpace(field, self.n+1)
        self.projection_points(V)
        
        # Make sure to trim off the identity element
        self.blocks = [self.projection([i for i in subspace][1:]) for subspace in V.subspaces(self.d+1)]

    def generate_hyperplanes(self):
        '''
        Generates the design isomorphic to the hyperplanes of the projective geometry
        '''
        # Uses Singer's method. Since there is an identification from an n+1 dimensional vector space over F_q to F_{q^(n+1)}, and
        # multiplication by x in F_{q^(n+1)} maps subspaces to subspaces, it is an automorphism of the design. Since we can identify
        # each element of the finite field Q(x) with a number such that x^i = Q(x), we can create an additive difference set to generate
        # our geometry with.

        field = GF(self.q**(self.n+1))
        #powers = field.primitive_element().powers(self.q**(self.n+1))
        powers = field.primitive_element().powers((self.q**(self.n+1)-1)/(self.q-1)-1)
        print(powers)
        p = field.base().order()
        alpha = log(self.q)/log(p)
        diffset = [i for i in range(len(powers)) if powers[i].polynomial().degree() <= self.n*alpha -1 ]
        print(diffset)
        diffdesign = difference_method([diffset], (self.q**(self.n+1)-1)/(self.q-1), check=False)
        diffdesign.generate()
        self.blocks = diffdesign.blocks
        self.V = diffdesign.V

class difference_method(BIBD):
    def __init__(self, difference_sets, n, repetition="False", check=True):
        self.difference_sets = difference_sets
        self.n = n
        self.repetition = repetition
        if check == True:
            self.verify()
        
    def verify(self):
        differences = [0 for i in range(self.n)]
        for diffset in self.difference_sets:
            for x,y in permutations(diffset, int(2)):  # No idea why int(2) is needed
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
        self.V = [designpoint(i) for i in range(self.n)]
        sets = []
        for difference_set in self.difference_sets:
            sets += [[self.V[(x+i)%self.n] for x in difference_set] for i in range(self.n)]
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

for v in range(5,20):
    for k in range(3,v):
        for lambduh in range(1,5):
            try:
                x = BIBD([v,k,lambduh])
                exist = x.existence()
                print(x.parameters, exist.exist, exist.message)
            except AssertionError:
                pass
