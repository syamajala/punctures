class Plat():

    def __init__ (self, minbraid):
        self.strands = minbraid.s*2;
        self.minbraid = minbraid;
        
    @classmethod
    def plat(cls, minbraid):
        return Plat(minbraid)

    def compute_plat():
        # we want to start from the right most strand
        # and move to the second strand
        i = self.strands
        j = 2
        ans = []
        

        while i != j:
            ans.append(chr(i+97-2))
            i = i - 1
                
        
