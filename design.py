from finitefield import *

def binom(n,r,q):
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
        return bino

class BIBD():
    def __init__(self, parameters):
        self.parameters = parameters
        self.v, self.k, self.lambduh = (int(parameters[0]), int(parameters[1]), int(parameters[2]))
    
        self.parameters = (self.v, self.k, self.lambduh)
        self.b = self.lambduh * self.v * (self.v-1)/(self.k*(self.k-1))
        self.r = self.lambduh * (self.v-1)/(self.k-1)
        assert(self.b.is_integer() and self.r.is_integer()), "inadmissible design"
        if self.b == self.v:
            self.symmetric = True
        else:
            self.symmetric = False

    def derived(self):
        if self.symmetric == True and self.lambduh > 1:
            return BIBD([self.k, self.lambduh, self.lambduh - 1])

    def residual(self):
        if self.symmetric == True:
            return BIBD([self.v-self.k, self.k - self.lambduh, self.lambduh])

    def complement(self):
        return BIBD([self.v, self.v-self.k, self.b-2*self.r+self.lambduh])

class AG(BIBD):
    def __init__(self,  n, q, d):
        self.n = n
        self.q = q
        self.d = d
        self.v, self.k, self.lambduh = (q**n,   q**d,  binom(n-1, d-1, q))
        self.parameters = [int(self.v), int(self.k), int(self.lambduh)]
        BIBD.__init__(self, self.parameters)

    def generate(self):
        field = GF(self.q)
        omega = field.gen()
        self.elements = [omega**i for i in range(self.q-1)]
    


class PG(BIBD):
    def __init__(self,  n, q, d):
        self.v, self.k, self.lambduh = (binom(n+1,  1, q), binom(n, 1, q), binom(n-1, 1, q))
        self.parameters = (int(self.v), int(self.k), int(self.lambduh))
        BIBD.__init__(self, self.parameters)

def genBIBD():
    for n in range(3,50):
        for k in range(1,n+1):
            for lambduh in range(1,5):
                try:
                    x = BIBD([n,k,lambduh])
                    print(x.parameters)
                    #print(x.parameters, x.complement().parameters, x.derived().parameters, x.residual().parameters)
                except AssertionError:
                    pass

def listAG():
    for n in range(2,5):
        for q in [2,3,4,5,8,9]:
            for d in range(1, n):
                try: 
                    x = AG(n,q,d)
                    print(n, " & ", q, " & ", d, " & ", x.parameters, " & ", x.complement().parameters, "\\\\")
                except AssertionError:
                    pass
                except AttributeError:
                    pass

def listPG():
    for n in range(2,5):
        for q in [2,3,4,5,8,9]:
            for d in range(1, n):
                try:
                    x = PG(n,q,d)
                    print(n, " & ", q, " & ", d, " & ", x.parameters, " & ", x.complement().parameters, " & ", x.derived().parameters, " & ", x.residual().parameters, "\\\\") 
                except AssertionError:
                    print(n, " & ", q, " & ", d, " & ", x.parameters, " & ", x.complement().parameters, " & N/A & ",  x.residual().parameters, "\\\\")
                except AttributeError:
                    print(n, " & ", q, " & ", d, " & ", x.parameters, " & ", x.complement().parameters, " & N/A & ", x.residual().parameters, "\\\\")

def listbinom():
    for q in [2,3,4,5,8,9]:
        for n in range(1,6):
            string = ""
            for r in range(0,n):
                string = str(r) + " & "
                string += str(int(binom(n,r,q))) + " & "
            for r in range(n,6):
                string += "- & "
            print(n, "&", string, "\\\\")

x = AG(5,11,3)
x.generate()
