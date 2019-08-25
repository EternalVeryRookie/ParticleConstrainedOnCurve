import numpy as np
from scipy.integrate import solve_ivp, Radau

from Spline2D import Spline2D

#x[0] = x, x[1] = dx/dt
#滑らかな2次元曲線に拘束された質点の運動方程式
def restraintedEOM(mass, g, resistance_coeff, restrainCurve):

  def equationOfMotion(t, x):
    C = restrainCurve.getValue(x[0])
    dCdx = restrainCurve.getDifferentialValue(x[0])
    d2Cdx2 = restrainCurve.get2OrderDifferentialValue(x[0])
    norm = dCdx[0]**2 + dCdx[1]**2
    sqrt_norm = np.sqrt(norm)
    speed = np.abs(x[1])
    dxdt = x[1]
    dvdt = (-mass*x[1]*x[1]*(dCdx[0]*d2Cdx2[0] + dCdx[1]*d2Cdx2[1]) + mass*g*dCdx[1] - resistance_coeff*(sqrt_norm**3)*speed*x[1]) / (mass*norm)

    return [dxdt, dvdt]

  return equationOfMotion


#曲線に拘束された質点のシミュレーター
#かかる力は重力と摩擦
#曲線はスプライン曲線
class MassPointRestraintedCurveSimulator:
  def __init__(self, mass=1):
    self.mass = mass
    self.g = 9.80665
    self.spline = Spline2D()


  def timeDevelop(self, time_range, init_conditions=[0,0]):
    start = time_range[0]
    end = time_range[len(time_range) - 1]

    s = solve_ivp(restraintedEOM(self.mass, self.g, 0.01, self.spline), [start, end], init_conditions, t_eval=time_range, method="RK45")
    return s.t, s.y
