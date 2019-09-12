from sympy import *

#3次曲線と点の距離を陽に書き下すプログラム

X = Symbol("X")
Y = Symbol("Y")

t = Symbol("t")

l_x = Symbol("l_x")
l_y = Symbol("l_y")

a_x = Symbol("a_x")
b_x = Symbol("b_x")
c_x = Symbol("c_x")
d_x = Symbol("d_x")

a_y = Symbol("a_y")
b_y = Symbol("b_y")
c_y = Symbol("c_y")
d_y = Symbol("d_y")

print( expand( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2 ) )
print( diff  ( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2, t) )
print("\n\n")

dist = (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2
diff = (-2*b_x - 2*c_x*(-2*l_x + 2*t) - 6*d_x*(-l_x + t)**2)*(X - a_x - b_x*(-l_x + t) - c_x*(-l_x + t)**2 - d_x*(-l_x + t)**3) + (-2*b_y - 2*c_y*(-2*l_y + 2*t) - 6*d_y*(-l_y + t)**2)*(Y - a_y - b_y*(-l_y + t) - c_y*(-l_y + t)**2 - d_y*(-l_y + t)**3)
a = Poly(expand(diff), t)
print(a.coeffs())
