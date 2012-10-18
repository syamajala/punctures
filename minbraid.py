class MinBraid():

   def __init__ (self, s, c):
      self.strands = s
      self.crossings = c
      self.braiduniverse = CartesianProduct(*[range(1, s)]*c)

   @classmethod
   def braid(cls, s, c):
      return MinBraid(s, c)

   def minbraid(self):
      self.filter1();
      self.fliter2();
      self.filter3();
      self.filter4();
      map(self.minbraid_help, self.braiduniverse)
      
   def minbraid_help(self, l):
      ans = []
      for i in l:
         if i%2 == 1:
            ans.append(chr(i+64))
         elif i%2 == 0:
            ans.append(chr(i+96))
      return ans

   def filter1(self):
      self.braiduniverse = [i for i in self.braiduniverse if i[0] == 1]

   def filter2(self):
      self.braiduniverse = filter(self.filter2_help, self.braiduniverse)

   def filter2_help(self, l):      
      s = range(1, self.strands)
      for i in s:
         if l.count(i) < 2:
            return false
         
      return true
   
   def filter3(self, l):
      pass

   def filter4(self, l):
      pass


      
