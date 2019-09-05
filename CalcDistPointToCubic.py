import sys
from sympy import *
#cubic : a + b*(t-l) + c*((t-l)**2) + d*((t-l)**3
#最近傍点、最近傍点までの距離、その点のパラメーターをリターンする
def calcDistPointToCubic(point, cubic):
  t = Symbol("t", real=True)

  X = point[0]
  Y = point[1]

  l_x = cubic[0].domain_lower
  l_y = cubic[1].domain_lower

  a_x = cubic[0].a
  b_x = cubic[0].b
  c_x = cubic[0].c
  d_x = cubic[0].d

  a_y = cubic[1].a
  b_y = cubic[1].b
  c_y = cubic[1].c
  d_y = cubic[1].d

  dist = (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2
  diff = (-2*b_x - 2*c_x*(-2*l_x + 2*t) - 6*d_x*(-l_x + t)**2)*(X - a_x - b_x*(-l_x + t) - c_x*(-l_x + t)**2 - d_x*(-l_x + t)**3) + (-2*b_y - 2*c_y*(-2*l_y + 2*t) - 6*d_y*(-l_y + t)**2)*(Y - a_y - b_y*(-l_y + t) - c_y*(-l_y + t)**2 - d_y*(-l_y + t)**3)
  sol = solve(diff)
  u_x = cubic[0].domain_upper
  min_dist = sys.float_info.max
  point = None
  solve_param = None
  for s in sol:
    if s <= u_x and s >= u_x:
      d = dist.subs(t, s)
      if min_dist > d:
        min_dist = d
        point = [cubic[0].get_value(s), cubic[1].get_value(s)]
        solve_param = s

  return [min_dist, point, solve_param]




#3次曲線と点の距離を陽に書き下すプログラム

X = Symbol("X")
Y = Symbol("Y")

t = Symbol("t", real=True)

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

print("距離の関数 : ", end="")
dist = expand( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2 )
print( expand( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2 ) )
print("上記関数の最小値を取る点が最近傍点。微分=0となる点を求め、それらを上記関数に代入し、実際に最小値を取るものを探す")

print("距離の関数の微分 : ", end="")
print( diff  ( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2, t) )


diff = diff  ( (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2, t)
#p = poly(diff)
#print(expand(diff))
#print(p.all_coeffs())

solve = solve(diff, t)
print(solve)

for s in [y for y in solve if (y >= l_x) and (y >= l_y)]:
  print(dist.subs([(t, s)]))
