import string

class MinBraid():

    def __init__ (self, s, c):
        self.strands = s
        self.crossings = c
        self.braiduniverse = CartesianProduct(*[range(1, s)]*c)

    @classmethod
    def braid(cls, s, c):
        return MinBraid(s, c)

    def minbraid(self):
        self.braiduniverse = self.filter1()
        self.braiduniverse = self.filter2()
        self.braiduniverse = self.filter3()
        self.braiduniverse = self.filter4()
        return self
        #        return map(self.minbraid_help, self.braiduniverse)
      
    def minbraid_help(self, l):
        ans = []
        for i in l:
            if i%2 == 1:
                ans.append(chr(i+64))
            elif i%2 == 0:
                ans.append(chr(i+96))
        return ans

    def filter1(self):
        return [i for i in self.braiduniverse if i[0] == 1]

    def filter2(self):
        return filter(self.filter2_help, self.braiduniverse)

    def filter2_help(self, l):      
        s = range(1, self.strands)
        for i in s:
            if l.count(i) < 2:
                return false
 
        return true
   
    def filter3(self):        
        return self.braiduniverse

    def filter4(self):
        return self.braiduniverse

    def smaller(self, b1, b2):
        if len(b1) != len(b2):
            return false

        for i in range(0, len(b1)):
            if b1[i] > b2[i]:
                return b2

        return b1

    def alternating(self, u):
        braid = self.minbraid_help(u)

        upper = [string.ascii_uppercase[i] for i in range(0, 26, 2)]
        lower = [string.ascii_lowercase[i] for i in range(1 ,26, 2)]

        for i in braid:
            if i in upper:
                if not self.odd(i):
                    return false
            elif i in lower:
                if not self.even(i):
                    return false
                
        return true
        
    def even(self, c):
        n = ord(string.lower(c))-97
        if n%2 == 0:
            return false
        else:
            return true

    def odd(self, c):
        return not self.even(c)
        
class Plat():

    def __init__ (self, minbraid):
        self.strands = minbraid.strands*2
        self.minbraid = minbraid
        
    @classmethod
    def plat(cls, minbraid):
        return Plat(minbraid)

    def compute_plat():
        
