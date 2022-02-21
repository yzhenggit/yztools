from george.modeling import Model
class Poly1D(Model):
    parameter_names = ("a0", "a1")

    def get_value(self, x):
        return self.a0 + self.a1*x
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        
class Poly2D(Model):
    parameter_names = ("a0", "a1", "a2")

    def get_value(self, x):
        return self.a0 + self.a1*x + self.a2*x**2
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        self.a2 = p[2]

class Poly3D(Model):
    parameter_names = ("a0", "a1", "a2", "a3")

    def get_value(self, x):
        return self.a0 + self.a1*x + \
               self.a2*x**2 + self.a3*x**3
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        self.a2 = p[2]
        self.a3 = p[3]
                
class Poly4D(Model):
    parameter_names = ("a0", "a1", "a2", "a3", "a4")

    def get_value(self, x):
        return self.a0 + self.a1*x + \
               self.a2*x**2 + self.a3*x**3 + \
               self.a4*x**4
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        self.a2 = p[2]
        self.a3 = p[3]
        self.a4 = p[4]
        
class Poly5D(Model):
    parameter_names = ("a0", "a1", "a2", "a3", "a4", "a5")

    def get_value(self, x):
        return self.a0 + self.a1*x + \
               self.a2*x**2 + self.a3*x**3 + \
               self.a4*x**4 + self.a5*x**5
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        self.a2 = p[2]
        self.a3 = p[3]
        self.a4 = p[4]
        self.a5 = p[5]

class Poly6D(Model):
    parameter_names = ("a0", "a1", "a2", "a3", "a4", "a5", "a6")

    def get_value(self, x):
        return self.a0 + self.a1*x + \
               self.a2*x**2 + self.a3*x**3 + \
               self.a4*x**4 + self.a5*x**5 + \
               self.a6*x**6
    def set_vector(self, p): 
        self.a0 = p[0]
        self.a1 = p[1]
        self.a2 = p[2]
        self.a3 = p[3]
        self.a4 = p[4]
        self.a5 = p[5]
        self.a6 = p[6]
