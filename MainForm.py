import tkinter as tk
from tkinter import *
import time
import numpy as np
import math
from copy import copy

from RestraintedEOM import MassPointRestraintedCurveSimulator

#canvas空間とシミュレーション空間を分けて考える
#canvas空間をそのままシミュレーションに利用すると扱う数値が大きくて誤差が大きくなるため
class MainForm(tk.Frame):
  def __init__(self, master=None, width=500, height=500):
    super().__init__(master)
    self.master = master
    self.pack()
    self.window_width = width + 200
    self.window_height = height + 50
    self.canvas_width = width
    self.canvas_height = height
    self.initWidgets()
    self.ctrl_p_radius = 10
    self.is_pick_ctrl_p = False
    self.pick_ctrl_p_index = -1
    self.is_simu_running = False
    self.max_ctrl_p_num = 10
    self.is_mouse_on_curve = False
    self.select_curve_index = -1
    self.dist_mouse_to_curve_th = 0.01#シミュレーション空間上での距離
    self.simulator = MassPointRestraintedCurveSimulator()
    self.addControlPoint([self.canvas_width - self.ctrl_p_radius, self.canvas_height - self.ctrl_p_radius])

    size = str(self.window_width)+"x"+str(self.window_height)
    self.master.geometry(size)

  def draw_canvas(self):
    self.canvas.delete("line")
    self.canvas.delete("ctrl_p")
    self.canvas.delete("mass_point")
    self.draw_curve()
    self.draw_ctrl_p()

  def draw_curve(self):
    points = self.simulator.spline.sampling(10)
    if self.is_mouse_on_curve:
      color = "green"
    else:
      color = "black"

    for i in range(len(points)-1):
      self.canvas.create_line(points[i][0]*self.canvas_width, points[i][1]*self.canvas_height, points[i+1][0]*self.canvas_width, points[i+1][1]*self.canvas_height, tag="line", fill=color, width=5)

  def draw_ctrl_p(self):
    ctrl_ps = self.simulator.spline.control_points
    color = "red"
    for p in reversed(ctrl_ps):
      p[0] *= self.canvas_width
      p[1] *= self.canvas_height
      self.canvas.create_oval(p[0]-self.ctrl_p_radius, p[1]-self.ctrl_p_radius, p[0]+self.ctrl_p_radius, p[1]+self.ctrl_p_radius, fill=color, tag="ctrl_p")
      color = "blue"

  def addControlPoint(self, point):
    if len(self.simulator.spline.control_points) - 1 >= self.max_ctrl_p_num:
      return

    point_copy = copy(point)
    point_copy[0] /= self.canvas_width
    point_copy[1] /= self.canvas_height
    self.simulator.spline.addControlPoint(point_copy)
    self.draw_canvas()

  #先頭に挿入する
  def insertControlPoint(self, point, index):
    if len(self.simulator.spline.control_points) - 1 >= self.max_ctrl_p_num:
      return

    point_copy = copy(point)
    point_copy[0] /= self.canvas_width
    point_copy[1] /= self.canvas_height
    self.simulator.spline.insertControlPoint(point_copy, index)
    self.draw_canvas()


  def pickCtrl(self, point_on_canvas):
    control_points = self.simulator.spline.control_points
    for index in range(len(control_points)-1): #最後の制御点は移動させない
      point = control_points[index]
      dx = point_on_canvas[0] - self.canvas_width*point[0]
      dy = point_on_canvas[1] - self.canvas_height*point[1]
      if (dx**2 + dy**2)< self.ctrl_p_radius**2:
        return index

    return -1

  def onLeftClick(self, evt):
    if self.is_simu_running:
      return

    self.pick_ctrl_p_index = self.pickCtrl([evt.x, evt.y])
    if self.pick_ctrl_p_index >= 0:
      return

    if self.is_mouse_on_curve:
      self.insertControlPoint([evt.x, evt.y], self.select_curve_index)
    else :
      self.insertControlPoint([evt.x, evt.y], 0)

  def onRightClick(self, evt):
    if self.is_simu_running:
      return

    control_points = self.simulator.spline.control_points
    for index in range(len(control_points)-1): #最後の制御点は消させない
      point = control_points[index]
      dx = evt.x - point[0] * self.canvas_width
      dy = evt.y - point[1] * self.canvas_height
      if (dx**2 + dy**2)< self.ctrl_p_radius**2:
        self.simulator.spline.removeControlPoint(index)
        self.draw_canvas()
        break


  def startSimulation(self):
    ctrl_ps = self.simulator.spline.control_points
    #制御点の数が2個未満のとき(曲線が生成されていないとき)は何もしない
    if len(ctrl_ps) < 2:
        return

    self.start_btn.config(state="disable")
    self.is_simu_running = True
    start_point = ctrl_ps[0]
    norm = self.simulator.spline.getDifferentialValue(0)
    norm = norm[0]**2 + norm[1]**2
    E = 0.01 #わずかに画面外に出れるようなエネルギーを与える
    U = -9.80665 * start_point[1]
    V = np.sqrt(2*(E-U)/norm)
    domain_of_def = [0, len(ctrl_ps) - 1]
    dt = 0.001
    init_condition = [0, V]
    #時間の単位は秒で統一
    sec_per_frame = 1/30
    elapsed_time = 0
    update_speed = 16
    while(True):
      if not(self.is_simu_running):
        break
      start_loop = time.perf_counter()
      for i in range(update_speed):
        _, solve = self.simulator.timeDevelop([elapsed_time, elapsed_time+dt], init_conditions=init_condition)
        s = solve[0][len(solve[0])-1]
        elapsed_time += dt

        #曲線の外に出ようとしたら座標を押し戻して速度を反転
        if s < domain_of_def[0]:
          s = domain_of_def[0]
          solve[1][len(solve[1])-1] *= -1
        elif s > domain_of_def[1]:
          s = domain_of_def[1]
          solve[1][len(solve[1])-1] *= -1
          break
        init_condition=[solve[0][len(solve[0])-1], solve[1][len(solve[1])-1]]

      p = self.simulator.spline.getValue(s)
      p[0] = int(p[0] * self.canvas_width)
      p[1] = int(p[1] * self.canvas_height)
      self.canvas.delete("mass_point")
      self.canvas.create_oval(p[0]-self.ctrl_p_radius, p[1]-self.ctrl_p_radius, p[0]+self.ctrl_p_radius, p[1]+self.ctrl_p_radius, fill="green", tag="mass_point")
      self.canvas.update()
      if time.perf_counter() - start_loop < sec_per_frame:
        time.sleep((sec_per_frame - (time.perf_counter() - start_loop))/1.1)
      else:
        print("処理落ち")
      self.elapsed_time_label["text"] = "{:.3f}".format(elapsed_time)
      if s == domain_of_def[1]:
        break

    self.is_simu_running = False
    self.start_btn.config(state="normal")


  def stopSimulation(self):
    self.is_simu_running = False


  def clearCtrlPs(self):
    if self.is_simu_running:
      return

    ctrl_ps = self.simulator.spline.control_points
    for i in range(len(ctrl_ps)-1):
      self.simulator.spline.removeControlPoint(0)
    self.draw_canvas()


  def onRelease(self, evt):
    self.pick_ctrl_p_index = -1

  def onDragg(self, evt):
    if self.is_simu_running:
      return

    if self.pick_ctrl_p_index < 0:
      return

    if evt.x < 0 or evt.y < 0 or evt.x > self.canvas_width or evt.y > self.canvas_height:
      return

    point = [evt.x/self.canvas_width, evt.y/self.canvas_height]
    self.simulator.spline.moveControlPoint(point, self.pick_ctrl_p_index)

    self.draw_canvas()


  def leave(self, evt):
    self.pick_ctrl_p_index = -1


  def mouseMove(self, evt):
    action = (lambda: 0)
    point = [evt.x/self.canvas_width, evt.y/self.canvas_height]
    d, point, param, min_dist_curve_index = self.simulator.spline.calcDistPointToSpline(point)
    th = 0.00001
    print(min_dist_curve_index)
    if d < th:
      if not(self.is_mouse_on_curve):
        action = self.draw_canvas
      self.is_mouse_on_curve = True
      self.select_curve_index = min_dist_curve_index + 1
    else:
      if self.is_mouse_on_curve:
        action = self.draw_curve
      self.is_mouse_on_curve = False

    action()


  def initWidgets(self):
    self.canvas = tk.Canvas(self, width=self.canvas_width, height=self.canvas_height, bd=2, bg="white")
    self.canvas.grid(column=0,row=0, rowspan=10)

    self.elapsed_time_label = tk.Label(self, text="0.000", width=10, font=("", 20))
    self.elapsed_time_label.grid(column=1, row=4)

    self.start_btn = tk.Button(self, text="スタート", bd=2, width=20, command=self.startSimulation)
    self.start_btn.grid(column=1, row=5)

    self.stop_btn = tk.Button(self, text="ストップ", bd=2, width=20, command=self.stopSimulation)
    self.stop_btn.grid(column=1, row=6)

    self.clear_ctrlps_btn = tk.Button(self, text="曲線クリア", bd=2, width=20, command=self.clearCtrlPs)
    self.clear_ctrlps_btn.grid(column=1, row=7)

    self.canvas.bind("<ButtonPress-1>", self.onLeftClick)
    self.canvas.bind("<ButtonPress-3>", self.onRightClick)
    self.canvas.bind("<ButtonRelease-1>", self.onRelease)
    self.canvas.bind("<B1-Motion>", self.onDragg)
    self.canvas.bind("<Motion>", self.mouseMove)
    self.canvas.bind("<Leave>", self.leave)





root = tk.Tk()
root.title("ParticleConstrainedOnCurve")
form = MainForm(root, 1200, 600)
form.mainloop()
