import sys
sys.path.append('/home/seshu/dev/BordProgSage/')

from bordered import *
from scoop import futures


class Pdisk(PMC):
    
    def __init__(self, punctures, matching=[], intervals=[], arcslid=false, slidto=0):
        self.punctures = punctures
        self.pts = 4*punctures-1
        self.intervals = []        
        self.arcslid = arcslid
        self.slidto = slidto
        
        if matching == []:
            PMC.__init__(self, punctures, self.compute_matching(punctures))
            x = 0
            y = 3
        
            for i in range(0, punctures):
                self.intervals.append(range(x, y))
                x = y
                y = y+3
                
            self.intervals.append(range(max(max(self.intervals))+1, self.pts+1))

        else:
            self.pmc = PMC.__init__(self, punctures, matching)
            self.intervals = intervals
        
        self.label = label(self)
                                                                                           
    def compute_matching(self, punctures):
        matching = []
        i = 0
        for x in range(0, (self.pts)-punctures, 3):
            matching.append((x, x+2))
        for x in range(1, (self.pts)-punctures, 3):
            matching.append((x, (self.pts)-i))
            i = i+1
        return matching
          
    def idempotents(self, spinc=0):
        """Returns a list of the idempotents for this pointed matched circle with g+i strands.
        Idempotents are given as lists of elements of matching."""
        idems = []
        for i in range(1, (self.genus*2)+1+spinc):
            idems.extend(Combinations(self.matching, i).list())
        return idems

    def mstrands(self, strands):
        ans = []
    
        for i,v in enumerate(strands):
            x, y = v
    
            if (y-x != 1) and (len(strands) == 1):
                return true    
            elif (y-x != 1) and (len(strands) != 1):
                return false
    
            if i == 0:
                ans.append(y)
            elif (v == strands[-1]):
                ans.append(x)    
            else:
                ans.append(x)
                ans.append(y)
        
        for i in ans:        
            if ans.count(i) != 2:
                return false

        return true

    def compute_mlen_idem(self, strands):    
        x, y = strands[0]
        l = y-x
        ans = []
        i = 1
        
        while l > 1:
            p = y-i
            q = self.matched_point(p)
            if p > q:
                ans.append((q, p))
            elif q > p:
                ans.append((p, q))
            l = l - 1
            i = i + 1

        return (strands[0], ans)

    def compute_basis(self, spinc=0):
        # i don't know if this is the most efficient way to do this, but it works...        
        stra = []
        strands = []
        basis = []
        mlen = []
        top = []        
                
        # we have two cases to deal with:
        # generating the basis for something after an arcslide
        # generating the basis for something not arcslid

        if __name__ == '__main__':

            if self.arcslid:
                
                # much of this code is the same as below, we just have 2 special intervals
                # to deal with. the top and the interval slid to, if the top == slidto there
                # is only one special interval to deal with and the code works the same as the non special case

                tstr = []
                toint = self.intervals[self.slidto]
                
                for i in self.intervals:
                    if (i != self.intervals[-1]) and (i != toint):
                        for j in Combinations(i, 2).list():
                            stra.append(tuple(j))
                    elif (i == self.intervals[-1]):        
                        for j in Combinations(i, 2).list():
                            top.append(tuple(j))
                    elif (i == toint):
                        for j in Combinations(i, 2).list():
                            tstr.append(tuple(j))

                for i in top:
                    x, y = i
                    if (y-x) == 1:
                        stra.append(i)  

                for i in tstr:
                    x, y = i
                    if (y-x) == 1:
                        stra.append(i)

                for i in Combinations(stra, 2).list():
                    if self.mstrands(i):
                        strands.append(i)

                for i in stra:
                    x, y = i
                    if (y-x) < 2:
                        strands.append([i])
                    elif (y-x) == 2:
                        mlen.append([i])

                for i in range(1, len(self.intervals[-1])):
                    for j in Combinations(top, i).list():
                        if self.mstrands(j) and (j not in mlen):
                            mlen.append(j)

                for i in range(1, len(toint)):
                    for j in Combinations(tstr, i).list():
                        if self.mstrands(j) and (j not in mlen):
                            mlen.append(j)     
                                               
            if not self.arcslid:
            
                # first we generate all the 1 moving strands possibilities
            
                for i in self.intervals:
                    if (i != self.intervals[-1]):
                        for j in Combinations(i, 2).list():
                            stra.append(tuple(j))
                    elif (i == self.intervals[-1]):
                        for j in Combinations(i, 2).list():
                            top.append(tuple(j))
                            
                for i in top:
                    x, y = i
                    if (y-x) == 1:
                        stra.append(i)                                            

                # then we generate the two 2 moving strand possibilities,
                # where the two moving strands are on consecutive points on an interval

                for i in Combinations(stra, 2).list():
                    if self.mstrands(i):
                        strands.append(i)

                # next we separate the 1 moving strands that are longer than 1
                # from the ones of length 1

                for i in stra:
                    x, y = i
                    if (y-x) < 2:
                        strands.append([i])
                    elif (y-x) >= 2:
                        mlen.append([i])

                # then we check the top for strands longer than 1
                # and for multiple strands on top of each other

                for i in range(1, self.punctures):
                    for j in Combinations(top, i).list():
                        if self.mstrands(j) and (j not in mlen):
                            mlen.append(j)

            for i in mlen:
                if len(i) != 1:
                    strands.append(i)

            mlen = [i for i in mlen if len(i) == 1]

            ids = self.idem_dict(spinc)

            args = [(self, i, ids, spinc) for i in strands]
            ans = list(futures.map(alg_element, args))
            

            for i in ans:
                basis.extend(i)

            mlen_ids = list(futures.map(self.compute_mlen_idem, mlen))
            mlen_ids = dict(mlen_ids)
            
            args = [(self, i, ids, spinc) for i in mlen]
            ans = list(futures.map(alg_element, args))

            for i in ans:
                for j in i:
                    if set(j.left_idem).intersection(set(mlen_ids[j.strands[0]])):
                        basis.append(j)                          
            
            args = [([], i) for i in self.idempotents()]
            ans = list(futures.map(self.pstrand_diagram, args))

            basis.extend(ans)


        self.basis[spinc]= list(set(basis))
        self.basis_computed[spinc] = True


    def pstrand_diagram(self, args):
        strands, id = args
        return Strand_Diagram(self, strands, left_idem=id)

    @classmethod
    def split_matching(cls, p):
        return Pdisk(p)
            
    def zero_type_D(self):
        "Returns the standard type D structure for an zero-framed handlebody of genus k."
                
        s = []
        gens = []
                     
        # build a list of generators using the (i, j) labeling of the arc diagram
        
        idx = [range(0, 2)]
        for i in range(0, (self.punctures/2)-1):
            idx.append(map(lambda x: x + 2, idx[-1]))                        
        
        for i in CartesianProduct(*idx):
            s.append(list(self.label.mpairs))
            for j in i:
                s[-1].append(self.label.nmpairs[j])
        
        # build a list of DGens out of the list of generators created above
                                    
        q = 0             
        names = map(chr, range(97, 123))   
        idems = map(self.build_idem, s)
        for i in idems:
            gens.append(DGen(names[q], self, i))
            q = q + 1                                 
                                                
        # build alpha for the differential

        a = []
        for i in range(1, len(self.intervals[-1])+1):
            if i%2 == 1:
                a.append([(self.label.label2[i], self.label.label2[i+1])])
            else:
                a.append([(self.label.label2[(i, 0)], self.label.label2[(i, 2)]), a[-1][0], (self.label.label2[(i-1, 0)], self.label.label2[(i-1, 2)])])
                
        # compute the differential for each generator, first check if we can 
        # multiply on terms of alpha on the right by any idempotents, then check the left
        
        sa = []
        ids = self.idem_dict()

        if __name__ == '__main__':
            
            a = [(self, i, ids, 0) for i in a]

            ans = list(futures.map(alg_element, a))
                                
            for i in ans:
                sa.extend(i)
            
            dgens = list(gens)
            args = [(i, dgens) for i in sa]
                    
            diffsr = list(futures.map(self.diffr, args))
            
            args = [(i, dgens) for i in diffsr if i != None]
            
            diffsl = list(futures.map(self.diffl, args))
            
        # build a dict with keys generators, and values lists of pairs (a, x)
        # there is probably a better way to do this, because this is where most of the 
        # time gets eaten up when doing computations
              
        diffs = {}

        for i in diffsl:
            x,y,z = i
            diffs.setdefault(x, []).extend([(y,z)])
                                                                            
        for i in gens:
            if i not in diffs:
                diffs[i] = []
        
        return PTypeDStr(self, gens, diffs)        
         
    def diffr(self, args):
        sa, gens = args
        for i in gens:
            if sa.left_idem==i.idem:
                return (i, sa)
            
    def diffl(self, args):
        sa, gens = args
        x, y = sa
        for i in gens:
            if y.right_idem==i.idem:
               return (x, y, i)                                                                                                                            
    
    def build_idem(self, gen):
        "Given a generator for BSD(Bn) return the corresponding idempotent"
        idem = []
        t = []
        
        for i in gen:
            for j in i:
                t.append(self.label.label2[j])
            t.sort()
            idem.append(tuple(t))
            t = []
        
        return idem                        
            
    def show(self):
        ints = []
        ys = -.5
        ye = 0
        
        for i in range(0, len(self.intervals)):
            ye = max(self.intervals[i])
            ints.append(line([(0, ys), (0, ye)], rgbcolor=(0,0,0)))
            ys = ye + .5        
        
        circ = ints.pop(0)
        pict = circ
        for i in ints:
            pict = pict + i
                
        for i in range(4*self.genus):
            pict = pict + line([(-.25,i),(.25,i)], rgbcolor=(.1,.1,.1))
        
        for pair in self.matching:
            j=self.matching.index(pair)
            offset = j/(2*self.genus)
            pict = pict+line([(-.3,pair[0]),(-.5-offset,pair[0]),(-.5-offset,pair[1]),(-.3,pair[1])], rgbcolor=(.3-j/(5*self.genus),.6+j/(5*self.genus),.3-j/(5*self.genus)))
        pict.show(axes=False, xmin=-2,xmax=2, ymin=-.5, ymax=4*self.genus-.5)
    
    def arcslide(self, i, j):

        ints = [list(x) for x in self.intervals]

        for interval in ints:
            if (i in interval) and (j in interval):
                pmc = PMC.arcslide(self, i, j)
                interval.remove(i)              
                break                                                   
            elif (interval == ints[-1]) and (i not in interval) and (j not in interval):
                raise Exception("error")
                
        # figure out what j is matched to              
        z = self.matched_point(j)
                                
        # find what interval z is in and what interval j is in        
        for r,s in enumerate(self.intervals):
            if z in s:
                w = r
            if j in s:
                v = r
                        
        # build the new intervals, we need to look at 4 cases
        # 1. underslides from any interval but the top
        # 2. overslides from any interval but the top
        # 3. underslides from the top
        # 4. overslides from the top
        
        for r,s in enumerate(ints):
            if v != (len(ints)-1):
               if self.is_underslide(i,j) and (r<w) and (r>v):
                   ints[r] = [x-1 for x in s]
               elif self.is_overslide(i,j) and (r<w) and (r>=v):
                   ints[r] = [x-1 for x in s]
               elif r == w:    
                   ints[r] = ([x-1 for x in s[:s.index(z)+1]] + s[s.index(z):])
            elif v == (len(ints)-1):
               if (r>w) and (r<v):
                   ints[r] = [x+1 for x in s]
               elif r==w:
                   ints[r] = (s[:s.index(z)+1] + [x+1 for x in s[s.index(z):]])
               elif self.is_underslide(i,j) and (r == v):
                   ints[r] = ([x+1 for x in s[:s.index(j)]] + s[s.index(j):])
               elif self.is_overslide(i,j) and (r == v):
                   ints[r] = ([x+1 for x in s[:s.index(j)+1]] + s[s.index(j)+1:])

        return Pdisk(self.punctures, pmc.matching, ints, arcslid=true, slidto=w)        
        
    def chords(self, spinc=0):
        "Return sum_(xi a chord)a(xi), in specified strands grading. Defaults to middle strands grading."
        arc_list = list()
        ids = self.idem_dict()        
        answer = self.zero()
        ans = []
        if __name__ == '__main__':
            
            for i in self.intervals:
                for j in CartesianProduct(i, i):
                    x, y = tuple(j)
                    if (x != y) and (x < y):
                        arc_list.append((self, [tuple(j)], ids, spinc))
            
            ans = list(futures.map(alg_element, arc_list))
            
        for i in ans:
            answer.extend(i)
            
        return answer            

    def idem_dict(self,spinc=0):    
        
        # build a dictionary of sets of idempotents so we can save time doing computations that involve idempotents
        
        a = dict()
        idems = self.idempotents(spinc)
        for i in idems:
            for j in i:
                if j not in a:
                    a[j] = [frozenset(i)]
                else:
                    a[j].append(frozenset(i))
        return a
        
    def opposite(self):
        "Returns the PMC obtained by orientation-reversing self."
        revmatch = list()
        ints = list(self.intervals)
        for (i, j) in self.matching:
            p = 7 - i
            q = 7 - j
            if p > q:
                revmatch.append((q, p))
            elif q > p:
                revmatch.append((p, q))
             
        for i, v in enumerate(ints):
            ints[i] = map(lambda x: 7-x, v)

        return Pdisk(self.genus, revmatch, ints) 
        
    def check_interval(self, x, y):        
        for i in self.intervals:
            if (x in i) and (y in i):
                return true                

        return false                        
                    
class label():
    
    def __init__(self, pdisk):
        self.label1 = {} # dict of {x: (i, j)}
        self.label2 = {} # dict of {(i, j): x}
        self.mpairs = [[]]
        self.nmpairs = []
        p = len(pdisk.intervals)-1
        q = 0
        
        # we want to go from a consecutive labeling of the points of an arc diagram, 
        # ie 0, 1, 2, 3, 4, 5, ..., to something like 0 => (2, 0), 1 => (2, 1), ...
        # we also want to build two dictionaries to easily be able to convert back and forth 
        # between the different labels
        
        for i in pdisk.intervals:
            if i != pdisk.intervals[-1]:
                for j in i:                  
                    if j%3 != 1:
                        self.mpairs[q].append((p, j%3))
                    elif j%3 == 1:
                        self.nmpairs.append([p, (p, j%3)])
                    if j%3 == 2:
                        self.mpairs.append([])
                        q = q + 1
                    self.label1[j] = (p, j%3)   
            p = p - 1
            
            if i == pdisk.intervals[-1]:
                p = 1
                for j in i:
                     self.label1[j] = p
                     p = p + 1 

        for i in self.label1:
            self.label2[self.label1[i]] = i
               
        self.mpairs.pop()
                                                             
def alg_element(args):
    "Returns list (sum) of all algebra elements in A(pmc) with the given strands, and any consistent idempotents."
    pmc, strands, idems, spinc = args  
    g = pmc.genus+spinc
    answer = list()
    if not isinstance(strands, list):
        raise Exception("Strands input to alg_element should be a list.")
        
    if not idems:            
        for i in pmc.idempotents(spinc):
            if l_idem_compat(pmc,strands,i):
                answer.append(Strand_Diagram(pmc,strands,left_idem=i))
    
    if not strands:
        return []      

    # the list of idempotents gets pretty long really fast. 
    # ex, for 4 punctures its 255 elements, 6 is 4095...
    # here is a trick to speed up computations, at the cost of memory   
    
    if idems:
        id = []
        for i in strands:
            x, y = i
            z = pmc.matched_point(x)
            if x < z:
                id.append(set(idems[(x, z)]))
            elif x > z:
                id.append(set(idems[(z, x)]))

        id = list(id.pop().intersection(*id))

        for i in id:
            if l_idem_compat(pmc, strands, list(i)):
                answer.append(Strand_Diagram(pmc,strands,left_idem=list(i)))
    
    if not answer:
        id = []
        for i in strands:
            x, y = i
            z = pmc.matched_point(x)
            if x < z:
                id.append((x,z))
            elif x > z:
                id.append((z,x))
                
        return AlgElt(Strand_Diagram(pmc, strands, left_idem=id))
                                        
    return AlgElt(answer, pmc) #Used to just return answer.
    
class PTypeDStr(TypeDStr):
    
    def mor_to_d(self, other):
        """Returns the chain complex of homomorphisms from self to other, where other is a type D structure.
        
        Examples:
        sage: infty_type_D(1).mor_to_d(infty_type_D(1))
        Digraph on 2 vertices
        sage: infty_type_D(2).mor_to_d(m_one_type_D(2))
        Digraph on 9 vertices
        sage: infty_type_D(2).mor_to_d(m_one_type_D(2)).homology()
        1
        """
        if self.pmc != other.pmc:
            raise Exception("Can only compute Mor's between type D structures over the same algebra.")

        alg_gens = self.pmc.alg_basis()

        if __name__== '__main__':
            
            args = []
            for i in CartesianProduct(range(0, len(self.basis)), range(0, len(alg_gens)), range(0, len(other.basis))).list():
                x,y,z = tuple(i)
                args.append((self.basis[x],alg_gens[y],other.basis[z]))

            gens = list(futures.map(self.generate_gens, args))

            ans = alg_gens
            mult = []
            for i in range(2, self.pmc.punctures + 2):                
                ans = list(futures.map(self.mult, ans))
                ans = reduce(list.__add__, ans, [])                        
                ans = list(set(ans))
                mult.extend(ans)
                                           
            args = []
            
            for i in CartesianProduct(range(0, len(self.basis)), range(0, len(mult)), range(0, len(other.basis))).list():
                x,y,z = tuple(i)
                args.append((self.basis[x],mult[y],other.basis[z]))

            ans = list(futures.map(self.generate_gens, args))

            gens.extend(ans)
                                                                                        
            while gens.count(None) != 0:
                gens.remove(None)
                                                                                                                   
            diffs = dict()            
            args = [(i, other, diffs) for i in gens]
#            pool.map_async(self.diff, args, chunksize=1)
            for i in args:
                self.diff(i)

            diffs = dict(diffs)

        return ChainCx(diffs)

    def mult(self, i):
        ans = []
        b = self.pmc.alg_basis()
        for j in b:
            r = i.r_multiply(j)
            l = i.l_multiply(j)
            if r and (r.strands != []) and (r not in b) and (r not in ans):
                ans.append(r)
            if l and (l.strands != []) and (l not in b) and (l not in ans):
                ans.append(l)
        return ans                

    def diff(self, args):
        f, other, diffs = args
        df = list()
        #d(f(x))
        #First differentiate y in (x,a,y)
        deltay = other.delta1(f[2])
        ady = AlgElt([f[1]])*deltay
        for (b,z) in ady.basis_expansion():
           df.append((f[0],b,z))
        #Now differentiate a.
        for b in f[1].differential():
           df.append((f[0],b,f[2]))
        #f(dx)
        for xp in self.basis:
            for (b,w) in self.delta1(xp).basis_expansion():
                 if w==f[0]:
                    if b*f[1]:
                        df.append((xp,b*f[1],f[2]))
        diffs[f] = reduce_mod_2(df)
        
    def generate_gens(self, args):
        x,y,z = args
        if (x.idem==y.left_idem) and (z.idem==y.right_idem):
            return args
            
class PUnderslide(Underslide):

    def complementary_idem(self, idem):
        "Given an idempotent of self.pmc_1, return the complementary idempotent of self.pmc_2."       
        x,y = idem
        
        ids = self.pmc_2.idempotents()
        ids = [i for i in ids if len(i) == (len(ids[-1])-1)]
        
        if (x!=self.b1) and (y!=self.b1):
            p = (self.r(x),self.r(y))                    
        elif (x==self.b1) or (y==self.b1): 
            if self.b1p > self.b2p:
                p = (self.b2p, self.b1p)
            elif self.b1p < self.b2p:
                p = (self.b1p, self.b2p)        

        for i in ids:
            if p not in i:
                return (idem,i)
    
    def build_comp_idem(self, arg):
        id, ids = arg
        ans = [ids[i] for i in id]                        
        return (id, list(ans.pop().intersection(*ans)))

    def almost_complementary_idem(self, args):
        i, mstrand, mstrandp, mstrando = args         
        x, y = i
        a = list(x)
        b = list(y)        
        if ((mstrando in a) or (mstrando in b)) and ((mstrand not in a) or (mstrandp not in b)):
            if () in a:
                a.remove(())
            if mstrando not in a:
                a.append(mstrando)
            if mstrando not in b:
                b.append(mstrando)
            if mstrand in a:            
                a.remove(mstrand)
            if mstrandp in b:
                b.remove(mstrandp)
            a.sort()
            b.sort()
            return (a, b)
                    
    def generate_gen_idems(self):
        """Generates a list of pairs of idempotents for generators for the DD module. Stored in self.gen_idems                                                                                                                               
        Currently restricts to the middle SpinC structure.                                                                                                                                                                                   
        """
        #Generators of type X (complementary)                                                                                                                                                                                                
        xidems = [([()], self.pmc_2.idempotents()[-1])]
        pmc1idems = [i[0] for i in self.pmc_1.idempotents() if len(i) == 1]        

        if __name__ == '__main__':

            
            # 
            # what we do is we build a dictionary of complementary idempotents for idempotents consisting of 1 term
            # then we walk through the list of the rest of the idempotents and for each element in that idempotent,
            # we look up what the complementary idempotent was in the dict, and then take the intersection of everything
            #
            
            ids = list(futures.map(self.complementary_idem, pmc1idems))

            for i in ids:
                x, y = i
                xidems.append(([x], y))
                                            
            ids = dict([(i,frozenset(v)) for i,v in ids])
                
            pmc1idems = [i for i in self.pmc_1.idempotents() if len(i) != 1]
            args = [(i, ids) for i in pmc1idems]
            
            xidems.extend(list(futures.map(self.build_comp_idem, args)))
                                            
        #Generators of type Y (sub-complementary)                                                                                                                                                                                            

            if self.b1 > self.b2:
                mstrand = (self.b2, self.b1)
            elif self.b1 < self.b2: # i don't think this can happen for an underslide, but better safe than sorry.
                mstrand = (self.b1, self.b2)
            if self.b1p > self.b2p:
                mstrandp = (self.b2p, self.b1p)
            elif self.b1p < self.b2p:
                mstrandp = (self.b1p, self.b2p)
            if self.c1 > self.c2:
                mstrando = (self.c2, self.c1)
            elif self.c2 > self.c1:
                mstrando = (self.c1, self.c2)    
                  
            args = [(i, mstrand, mstrandp, mstrando) for i in xidems]      

            yidems = futures.map(self.almost_complementary_idem, args)
            
            while yidems.count(None) != 0:
                yidems.remove(None)
             
            for i in yidems:
                while yidems.count(i) != 1:
                    yidems.remove(i)
                                                                                                            
            self.gen_idems = xidems + yidems
            
        
    def dd_mod(self):
        "Return the type DD module for this arcslide."
        if self.ddstr != None:
            return self.ddstr
        if self.generators == None:
            self.generate_generators()
        self.generate_chords()
        diffs = dict()
        for x in self.generators:
            dx = DDElt({}, self.pmc_1, self.pmc_2rev)
            for y in self.generators:
                for (a,b) in self.allchords:
                    #If a is compatible with x's left idempotents and y's left idempotents, and b is compatible with x's right idempotents and y's right idempotents, add corresponding elt to dx.
                    if l_idem_compat(self.pmc_1, a, x.idem_1) and l_idem_compat(self.pmc_2rev, self.opposite_strands(b), x.idem_2):
                        aalgelt = AlgElt(Strand_Diagram(self.pmc_1, a, x.idem_1))*AlgElt(Strand_Diagram(self.pmc_1, [], y.idem_1))
                        balgelt = AlgElt(Strand_Diagram(self.pmc_2rev, self.opposite_strands(b), x.idem_2))*AlgElt(Strand_Diagram(self.pmc_2rev, [], y.idem_2))
                        abalgelt = aalgelt**balgelt
                        if abalgelt:
                            dx = dx + DDElt({y:abalgelt}, self.pmc_1, self.pmc_2rev)
            diffs[x] = dx
        self.ddstr = PTypeDDStr(self.pmc_1, self.pmc_2rev, self.generators, diffs)
        return self.ddstr
        
class POverslide(Overslide):
    
    def complementary_idem(self, idem):
        "Given an idempotent of self.pmc_1, return the complementary idempotent of self.pmc_2."       
        x,y = idem
        
        ids = self.pmc_2.idempotents()
        ids = [i for i in ids if len(i) == (len(ids[-1])-1)]
        
        if (x!=self.b1) and (y!=self.b1):
            p = (self.r(x),self.r(y))                    
        elif (x==self.b1) or (y==self.b1): 
            if self.b1p > self.b2p:
                p = (self.b2p, self.b1p)
            elif self.b1p < self.b2p:
                p = (self.b1p, self.b2p)        

        for i in ids:
            if p not in i:
                return (idem,i) 

    def build_comp_idem(self, arg):
        id, ids = arg
        ans = [ids[i] for i in id]                       
        return (id, list(ans.pop().intersection(*ans)))

    def almost_complementary_idem(self, args):
        i, mstrand, mstrandp, mstrando, mstrandop = args         
        x, y = i
        a = list(x)
        b = list(y)        
        if ((mstrando in a) or (mstrandop in b)) and ((mstrand not in a) or (mstrandp not in b)):
            if () in a:
                a.remove(())
            if mstrando not in a:
                a.append(mstrando)
            if mstrandop not in b:
                b.append(mstrandop)
            if mstrand in a:            
                a.remove(mstrand)
            if mstrandp in b:
                b.remove(mstrandp)
            a.sort()
            b.sort()
            return (a, b)
                    
    def generate_gen_idems(self):
        """Generates a list of pairs of idempotents for generators for the DD module. Stored in self.gen_idems                                                                                                                               
        Currently restricts to the middle SpinC structure.                                                                                                                                                                                   
        """
        #Generators of type X (complementary)                                                                                                                                                                                                
        xidems = [([()], self.pmc_2.idempotents()[-1])]
        pmc1idems = [i[0] for i in self.pmc_1.idempotents() if len(i) == 1]        

        if __name__ == '__main__':
            # 
            # what we do is we build a dictionary of complementary idempotents for idempotents consisting of 1 term
            # then we walk through the list of the rest of the idempotents and for each element in that idempotent,
            # we look up what the complementary idempotent was in the dict, and then take the intersection of everything
            #
            
            ids = list(futures.map(self.complementary_idem, pmc1idems))

            for i in ids:
                x, y = i
                xidems.append(([x], y))
                                            
            ids = dict([(i,frozenset(v)) for i,v in ids])
                
            pmc1idems = [i for i in self.pmc_1.idempotents() if len(i) != 1]
            args = [(i, ids) for i in pmc1idems]
            
            xidems.extend(list(futures.map(self.build_comp_idem, args)))
                                            
        #Generators of type Y (sub-complementary)                                                                                                                                                                                            

            # copied from PUnderslide.generate_gen_idems()
            
            if self.b1 > self.b2:
                mstrand = (self.b2, self.b1)
            elif self.b1 < self.b2: # i don't think this can happen for an underslide, but better safe than sorry.
                mstrand = (self.b1, self.b2)
            if self.b1p > self.b2p:
                mstrandp = (self.b2p, self.b1p)
            elif self.b1p < self.b2p:
                mstrandp = (self.b1p, self.b2p)
            if self.c1 > self.c2:
                mstrando = (self.c2, self.c1)
            elif self.c2 > self.c1:
                mstrando = (self.c1, self.c2)    
            if self.c1p > self.c2p:
                mstrandop = (self.c2p, self.c1p)
            elif self.c2p > self.c1p:
                mstrandop = (self.c1p, self.c2p)
                  
            args = [(i, mstrand, mstrandp, mstrando, mstrandop) for i in xidems]      

            yidems = list(futures.map(self.almost_complementary_idem, args))
            
            while yidems.count(None) != 0:
                yidems.remove(None)
            
            for i in yidems:
                while yidems.count(i) != 1:
                    yidems.remove(i)
                                                                            
            self.gen_idems = xidems + yidems
            

    def dd_mod(self):
        "Return the type DD module for this arcslide."
        if self.ddstr != None:
            return self.ddstr
        if self.generators == None:
            self.generate_generators()
        self.generate_chords()
        self.basic_choice()
        diffs = dict()
        for x in self.generators:
            dx = DDElt({}, self.pmc_1, self.pmc_2rev)
            for y in self.generators:
                #Process the determinate chords.
                for (a,b) in self.DetChords+self.choice_chords:
                    #If a is compatible with x's left idempotents and y's left idempotents, and b is compatible with x's right idempotents and y's right idempotents, add corresponding elt to dx.
                    if l_idem_compat(self.pmc_1, a, x.idem_1) and l_idem_compat(self.pmc_2rev, self.opposite_strands(b), x.idem_2):
                        aalgelt = AlgElt(Strand_Diagram(self.pmc_1, a, x.idem_1))*AlgElt(Strand_Diagram(self.pmc_1, [], y.idem_1))
                        balgelt = AlgElt(Strand_Diagram(self.pmc_2rev, self.opposite_strands(b), x.idem_2))*AlgElt(Strand_Diagram(self.pmc_2rev, [], y.idem_2))
                        abalgelt = aalgelt**balgelt
                        if abalgelt:
                            dx = dx + DDElt({y:abalgelt}, self.pmc_1, self.pmc_2rev)
                #Now process the indeterminate chords -- write.
            diffs[x] = dx
        self.ddstr = PTypeDDStr(self.pmc_1, self.pmc_2rev, self.generators, diffs)
        return self.ddstr            
            
class PArcslide(Arcslide):
    
    def __init__(self, pmc, i, j=i+1):
        self.pmc_1 = pmc
        self.pmc_2 = pmc.arcslide(i,j)
        self.i=i
        self.j=j
        if self.pmc_1.is_overslide(i,j):
            self.arcslide = POverslide(self.pmc_1,i,j)
        if self.pmc_1.is_underslide(i,j):
            self.arcslide = PUnderslide(self.pmc_1,i,j)
        self.pmc_2rev = self.arcslide.pmc_2rev
        self.b1 = self.arcslide.b1
        self.b1p = self.arcslide.b1p
        self.c1 = self.arcslide.c1
        self.c1p = self.arcslide.c1p
        self.c2 = self.arcslide.c2
        self.c2p = self.arcslide.c2p

class PTypeDDStr(TypeDDStr):
    
    def mor_to_d(self, other):
        """Returns the chain complex of homomorphisms over A(pmc_1) from self to other, where other is a type D structure.

        Examples:
        """
        if self.pmc_1 != other.pmc:
            raise Exception("Can only compute Mor's between structures over the same algebra.")
        #Compile a basis for Mor
        mor_basis = list()
        args = []
        
        for i in CartesianProduct(range(0, len(self.basis)), range(0, len(other.basis)), range(0, len(self.pmc_1.alg_basis()))).list():
            m, n, a = tuple(i)
            args.append((self.basis[m], other.basis[n], self.pmc_1.alg_basis()[a]))
        
        if __name__ == '__main__':

            
            mor_basis = list(futures.map(self.compute_basis, args))

            ans = self.pmc_1.alg_basis()
            mult = []
            for i in range(2, self.pmc_1.punctures + 2):                
                ans = list(futures.map(self.mult, ans))
                ans = reduce(list.__add__, ans, [])                        
                ans = list(set(ans))
                mult.extend(ans)
                                           
            args = []            
            for i in CartesianProduct(range(0, len(self.basis)), range(0, len(other.basis)), range(0, len(mult))).list():
                m, n, a = tuple(i)
                args.append((self.basis[m], other.basis[n], mult[a]))

            mor_basis.extend(list(futures.map(self.compute_basis, args)))
            
            while mor_basis.count(None) != 0:
                mor_basis.remove(None)

            diffs = dict()
            #delta^1 on Mor in three parts.
            #First: differentiating the algebra element
            for f in mor_basis:
                diffs[f]=DElt({},pmc=self.pmc_2)
                (m,a,n)=f.name
                for b in a.differential():
                    diffs[f]=diffs[f]+DElt({DGen((m,b,n),pmc=self.pmc_2,idem=m.idem_2):AlgElt([Strand_Diagram(self.pmc_2,[],m.idem_2)])})
            
            #Second: differentiating the image (element of other)
            #        print diffs   
            for f in mor_basis:   
                (m,a,n)=f.name
                adn = AlgElt([a])*other.delta1(n)
                for (b,p) in adn.basis_expansion():
                    diffs[f] = diffs[f]+DElt({DGen((m,b,p),pmc=self.pmc_2,idem=m.idem_2):AlgElt([Strand_Diagram(self.pmc_2,[],m.idem_2)])})

            #Third: differentiating the source (element of self). This is the only part which outputs non-idempotent algebra elements.
            #        print diffs
            for f in mor_basis:
                (m,a,n)=f.name
                for l in self.basis:
                    dl=self.delta1(l)
                    if m in dl.keys():
                        for b in dl[m].keys():
                            if b*a:
                                diffs[f] = diffs[f]+DElt({DGen((l,b*a,n),pmc=self.pmc_2,idem=l.idem_2):dl[m][b]})

            #        print diffs
            #Now view PMC as over opposite algebra:
            #First, reverse the idempotents.
            rev_basis = list()
            for f in mor_basis:
                (m,a,n)=f.name
                rev_idem = list()
                for (i,j) in m.idem_2:
                    rev_idem.append((4*self.pmc_2.genus-1-j,4*self.pmc_2.genus-1-i))
                rev_basis.append(DGen((m,a,n),pmc=self.pmc_2,idem=rev_idem))
            rev_diffs = dict()
            #Now, reverse the algebra element outputs.
            for f in diffs.keys():
                rev_diffs[f]=diffs[f].opposite()
        return PTypeDStr(self.pmc_2.opposite(), rev_basis, rev_diffs)        

    def compute_basis(self, args):
        m, n, a = args
        if (Set(m.idem_1)==Set(a.left_idem)) and (Set(n.idem)==Set(a.right_idem)):
            return DGen((m,a,n),pmc=self.pmc_2,idem=m.idem_2)
            
    def mult(self, i):
        ans = []
        b = self.pmc_1.alg_basis()
        for j in b:
            r = i.r_multiply(j)
            l = i.l_multiply(j)
            if r and (r.strands != []) and (r not in b) and (r not in ans):
                ans.append(r)
            if l and (l.strands != []) and (l not in b) and (l not in ans):
                ans.append(l)
        return ans
