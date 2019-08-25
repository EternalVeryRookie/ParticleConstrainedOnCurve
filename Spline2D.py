import numpy as np
import math
from copy import copy, deepcopy

class CubicFunction :
  def __init__(self, a, b, c, d, lower, upper):
    self.a = a
    self.b = b
    self.c = c
    self.d = d
    self.domain_lower = lower
    self.domain_upper = upper

  def get_value(self, t):
    """
    if  t < self.domain_lower or t > self.domain_upper:
      print("caution : Out of definition :" + str(t) + "\n")
    """

    variable = t - self.domain_lower
    return self.a + self.b * variable + self.c * pow(variable, 2) + self.d * pow(variable, 3)

  def get_differential_value(self, t):
    """
    if t < self.domain_lower or t > self.domain_upper:
      print("caution : Out of definition :" + str(t) + "\n")
    """

    variable = t - self.domain_lower
    return self.b + 2*self.c*variable + 3*self.d*variable*variable

  def get_2_order_differential_value(self, t):
    """
    if t < self.domain_lower or t > self.domain_upper:
      print("caution : Out of definition :" + str(t) + "\n")
    """

    variable = t - self.domain_lower
    return 2*self.c + 6*self.d*variable


#制御点間の媒介変数の幅は1
#制御点がN個あれば媒介変数の定義域は0<=t<=N-1
class Spline2D:
  def __init__(self):
    self.__control_points = []
    self.__curve_x = []
    self.__curve_y = []

  def addControlPoint(self, point):
    self.__control_points.append(copy(point))
    if len(self.__control_points) == 1:
      return
    lower = len(self.__control_points) - 2
    upper = lower + 1
    self.__curve_x.append(CubicFunction(0, 0, 0, 0, lower, upper))
    self.__curve_y.append(CubicFunction(0, 0, 0, 0, lower, upper))
    self.__calcSplineCoefficient()

  def insertControlPoint(self, point, index):
    if len(self.__control_points) == 0 or index < 0 or len(self.__control_points) <= index:
      return

    self.__control_points.insert(index, copy(point))
    lower = len(self.__control_points) - 2
    upper = lower + 1
    self.__curve_x.append(CubicFunction(0, 0, 0, 0, lower, upper))
    self.__curve_y.append(CubicFunction(0, 0, 0, 0, lower, upper))
    self.__calcSplineCoefficient()

  def removeControlPoint(self, index):
    if len(self.__control_points) <= 0 or len(self.__control_points) <= index:
      return

    self.__control_points.pop(index)
    if len(self.__curve_x) > 0:
      self.__curve_x.pop(len(self.__curve_x)-1)
      self.__curve_y.pop(len(self.__curve_y)-1)
      for i in range(len(self.__curve_x)): #定義域を再設定
        self.__curve_x[i].domain_lower = i
        self.__curve_x[i].domain_upper = i + 1

    self.__calcSplineCoefficient()

  def moveControlPoint(self, coordinate, index):
    if len(self.__control_points) <= 0 or len(self.__control_points) <= index:
      return

    self.__control_points[index][0] = coordinate[0]
    self.__control_points[index][1] = coordinate[1]
    self.__calcSplineCoefficient()

  def __calcSplineCoefficient(self):
    if len(self.__control_points) < 2:
      return

    curves_xy = [self.__curve_x, self.__curve_y]
    left_side_matrix = self.__calcLeftSideMatrix()
    right_side_vectors = self.__calcRightSideVectors()
    for dim in range(2):
      solve = np.linalg.solve(left_side_matrix, right_side_vectors[dim])

      for curveI in range(len(curves_xy[dim])):
        curves_xy[dim][curveI].a = solve[4 * curveI + 0]
        curves_xy[dim][curveI].b = solve[4 * curveI + 1]
        curves_xy[dim][curveI].c = solve[4 * curveI + 2]
        curves_xy[dim][curveI].d = solve[4 * curveI + 3]



  def __calcLeftSideMatrix(self):
    mat = np.zeros(( 4*len(self.__curve_x), 4*len(self.__curve_x) ))


    domainWidth = 1.0 #この実装では各3次曲線の定義域の幅は常に1
    for curveI in range(len(self.__curve_x)-1):

      #start point
      mat[4 * curveI][4 * curveI] = 1.0

      #end point
      for i in range(4):
      	mat[4 * curveI + 1][4 * curveI + i] = pow(domainWidth, i)


      mat[4 * curveI + 2][4 * curveI + 1] =  1.0
      mat[4 * curveI + 2][4 * curveI + 2] =  2.0 * (domainWidth)
      mat[4 * curveI + 2][4 * curveI + 3] =  3.0 * pow(domainWidth, 2)
      mat[4 * curveI + 2][4 * curveI + 5] = -1.0


      mat[4 * curveI + 3][4 * curveI + 2] = 2.0
      mat[4 * curveI + 3][4 * curveI + 3] = 6.0 * domainWidth
      mat[4 * curveI + 3][4 * curveI + 6] = -2.0


    curveN = len(self.__curve_x)-1
    #c_0 = 0
    mat[4 * curveN][2] = 1.0

    #start point
    mat[4 * curveN + 1][4 * curveN] = 1.0

    #end point
    mat[4 * curveN + 2][4 * curveN] = 1.0
    mat[4 * curveN + 2][4 * curveN + 1] = domainWidth
    mat[4 * curveN + 2][4 * curveN + 2] = pow(domainWidth, 2)
    mat[4 * curveN + 2][4 * curveN + 3] = pow(domainWidth, 3)

    #2������ = 0
    mat[4 * curveN + 3][4 * curveN + 2] = 2.0
    mat[4 * curveN + 3][4 * curveN + 3] = 6.0 * domainWidth

    return mat


  def __calcRightSideVectors(self):
    v_xy = np.zeros(( 2, 4*len(self.__curve_x) ))

    for dim in  range(2):
      for curveI in range(len(self.__curve_x)-1):
        v_xy[dim][4 * curveI + 0] = self.__control_points[curveI][dim] #curveIのstartPoint
        v_xy[dim][4 * curveI + 1] = self.__control_points[curveI+1][dim] #curveIのendPoint
        v_xy[dim][4 * curveI + 2] = v_xy[dim][4 * curveI + 3] = 0.0

      curveN = len(self.__curve_x) - 1
      #c_0 = 0
      v_xy[dim][4 * curveN] = 0

      v_xy[dim][4 * curveN + 1] = self.__control_points[curveN][dim]
      v_xy[dim][4 * curveN + 2] = self.__control_points[curveN+1][dim]

      v_xy[dim][4 * curveN + 3] = 0

    return v_xy

  def sampling(self, num_onesection):
    if len(self.__control_points) < 2:
      return []

    points = [None] * (num_onesection * len(self.__curve_x) + 1)
    for curveI in range(len(self.__curve_x)):
      for i in range(num_onesection):
        x = self.__curve_x[curveI].get_value(curveI + i/num_onesection)
        y = self.__curve_y[curveI].get_value(curveI + i/num_onesection)
        points[i + curveI * num_onesection] = [x, y]

    x = self.__curve_x[len(self.__curve_x)-1].get_value(len(self.__curve_x))
    y = self.__curve_y[len(self.__curve_y)-1].get_value(len(self.__curve_y))
    points[len(points)-1] = [x, y]
    return points

  def getValue(self, t):
    if len(self.__curve_x) == 0:
      return [None, None]

    if t < 0:
      index = 0
    elif t >= len(self.__curve_x):
      index = len(self.__curve_x) - 1
    else:
      index = math.floor(t)

    return [self.__curve_x[index].get_value(t), self.__curve_y[index].get_value(t)]

  def getDifferentialValue(self, t):
    if len(self.__curve_x) == 0:
      return [None, None]

    if t < 0:
      index = 0
    elif t >= len(self.__curve_x):
      index = len(self.__curve_x) - 1
    else:
      index = math.floor(t)

    return [self.__curve_x[index].get_differential_value(t), self.__curve_y[index].get_differential_value(t)]

  def get2OrderDifferentialValue(self, t):
    if len(self.__curve_x) == 0:
      return [None, None]

    if t < 0:
      index = 0
    elif t >= len(self.__curve_x):
      index = len(self.__curve_x) - 1
    else:
      index = math.floor(t)

    return [self.__curve_x[index].get_2_order_differential_value(t), self.__curve_y[index].get_2_order_differential_value(t)]




  @property
  def control_points(self):
    return deepcopy(self.__control_points)

if __name__ == "__main__":
  curve = Spline2D()
  curve.addControlPoint([1,2])
  curve.addControlPoint([1,2])
  curve.addControlPoint([1,2])
  curve.addControlPoint([1,2])
