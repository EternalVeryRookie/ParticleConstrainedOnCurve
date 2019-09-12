import sys
import numpy as np
#cubic : a + b*(t-l) + c*((t-l)**2) + d*((t-l)**3
#最近傍点、最近傍点までの距離、その点のパラメーターをリターンする
def calcDistPointToCubic(point, cubic):
  #t = Symbol("t", real=True)

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

  dist = lambda t: (X - (a_x + b_x*(t-l_x) + c_x*((t-l_x)**2) + d_x*((t-l_x)**3) ) )**2 + (Y - (a_y + b_y*(t-l_y) + c_y*((t-l_y)**2) + d_y*((t-l_y)**3)) )**2
  #diff = lambda t: (-2*b_x - 2*c_x*(-2*l_x + 2*t) - 6*d_x*(-l_x + t)**2)*(X - a_x - b_x*(-l_x + t) - c_x*(-l_x + t)**2 - d_x*(-l_x + t)**3) + (-2*b_y - 2*c_y*(-2*l_y + 2*t) - 6*d_y*(-l_y + t)**2)*(Y - a_y - b_y*(-l_y + t) - c_y*(-l_y + t)**2 - d_y*(-l_y + t)**3)

  #曲線との距離の式を微分した式を展開した時の各次数における係数の配列
  #事前にsympyを利用して求めたものをコピペ(リアルタイムに計算させると遅い)
  coeff = [6*d_x**2 + 6*d_y**2, 10*c_x*d_x + 10*c_y*d_y - 30*d_x**2*l_x - 30*d_y**2*l_y, 8*b_x*d_x + 8*b_y*d_y + 4*c_x**2 - 40*c_x*d_x*l_x + 4*c_y**2 - 40*c_y*d_y*l_y + 60*d_x**2*l_x**2 + 60*d_y**2*l_y**2, -6*X*d_x - 6*Y*d_y + 6*a_x*d_x + 6*a_y*d_y + 6*b_x*c_x - 24*b_x*d_x*l_x + 6*b_y*c_y - 24*b_y*d_y*l_y - 12*c_x**2*l_x + 60*c_x*d_x*l_x**2 - 12*c_y**2*l_y + 60*c_y*d_y*l_y**2 - 60*d_x**2*l_x**3 - 60*d_y**2*l_y**3, -4*X*c_x + 12*X*d_x*l_x - 4*Y*c_y + 12*Y*d_y*l_y + 4*a_x*c_x - 12*a_x*d_x*l_x + 4*a_y*c_y - 12*a_y*d_y*l_y + 2*b_x**2 - 12*b_x*c_x*l_x + 24*b_x*d_x*l_x**2 + 2*b_y**2 - 12*b_y*c_y*l_y + 24*b_y*d_y*l_y**2 + 12*c_x**2*l_x**2 - 40*c_x*d_x*l_x**3 + 12*c_y**2*l_y**2 - 40*c_y*d_y*l_y**3 + 30*d_x**2*l_x**4 + 30*d_y**2*l_y**4, -2*X*b_x + 4*X*c_x*l_x - 6*X*d_x*l_x**2 - 2*Y*b_y + 4*Y*c_y*l_y - 6*Y*d_y*l_y**2 + 2*a_x*b_x - 4*a_x*c_x*l_x + 6*a_x*d_x*l_x**2 + 2*a_y*b_y - 4*a_y*c_y*l_y + 6*a_y*d_y*l_y**2 - 2*b_x**2*l_x + 6*b_x*c_x*l_x**2 - 8*b_x*d_x*l_x**3 - 2*b_y**2*l_y + 6*b_y*c_y*l_y**2 - 8*b_y*d_y*l_y**3 - 4*c_x**2*l_x**3 + 10*c_x*d_x*l_x**4 - 4*c_y**2*l_y**3 + 10*c_y*d_y*l_y**4 - 6*d_x**2*l_x**5 - 6*d_y**2*l_y**5]

  sol = np.poly1d(coeff).roots
  u_x = cubic[0].domain_upper
  l_x = cubic[0].domain_lower
  min_dist = sys.float_info.max
  point = None
  solve_param = None
  for s in sol:

    if np.isreal(s) and s <= u_x and s >= l_x:
      d = dist(s)
      if min_dist > d:
        min_dist = d
        point = [cubic[0].get_value(s), cubic[1].get_value(s)]
        solve_param = s

  return [min_dist, point, solve_param]
