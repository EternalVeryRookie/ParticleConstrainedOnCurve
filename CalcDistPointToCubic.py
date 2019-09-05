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
