#!/usr/bin/env python
import numpy as np

def normalize(v):
  norm = np.linalg.norm(v)
  assert norm > 0
  return v / norm

def norm(v):
  return np.linalg.norm(v)

class Polynomial:
  def __init__(self, p):
    self.p = p

  def stretchtime(self, factor):
    recip = 1.0 / factor;
    scale = recip
    for i in range(1, len(self.p)):
      self.p[i] *= scale
      scale *= recip

  # evaluate a polynomial using horner's rule
  def eval(self, t):
    assert t >= 0
    x = 0.0
    for i in range(0, len(self.p)):
      x = x * t + self.p[len(self.p) - 1 - i]
    return x

  # compute and return derivative
  def derivative(self):
    return Polynomial([(i+1) * self.p[i+1] for i in range(0, len(self.p) - 1)])


class TrajectoryOutput:
  def __init__(self):
    self.pos = None   # position [m]
    self.vel = None   # velocity [m/s]
    self.acc = None   # acceleration [m/s^2]
    self.omega = None # angular velocity [rad/s]
    self.yaw = None   # yaw angle [rad]
    self.roll = None  # required roll angle [rad]
    self.pitch = None # required pitch angle [rad]
    self.jerk = None
    self.snap = None
# 4d single polynomial piece for x-y-z-yaw, includes duration.
class Polynomial4D:
  def __init__(self, duration, px, py, pz, pyaw):
    self.duration = duration
    self.px = Polynomial(px)
    self.py = Polynomial(py)
    self.pz = Polynomial(pz)
    self.pyaw = Polynomial(pyaw)

  # compute and return derivative
  def derivative(self):
    return Polynomial4D(
      self.duration,
      self.px.derivative().p,
      self.py.derivative().p,
      self.pz.derivative().p,
      self.pyaw.derivative().p)

  def stretchtime(self, factor):
    self.duration *= factor
    self.px.stretchtime(factor)
    self.py.stretchtime(factor)
    self.pz.stretchtime(factor)
    self.pyaw.stretchtime(factor)

  # see Daniel Mellinger, Vijay Kumar:
  #     Minimum snap trajectory generation and control for quadrotors. ICRA 2011: 2520-2525
  #     section III. DIFFERENTIAL FLATNESS
  def eval(self, t):
    result = TrajectoryOutput()
    # flat variables
    result.pos = np.array([self.px.eval(t), self.py.eval(t), self.pz.eval(t)])
    result.yaw = self.pyaw.eval(t)

    # 1st derivative
    derivative = self.derivative()
    result.vel = np.array([derivative.px.eval(t), derivative.py.eval(t), derivative.pz.eval(t)])
    dyaw = derivative.pyaw.eval(t)

    # 2nd derivative
    derivative2 = derivative.derivative()
    result.acc = np.array([derivative2.px.eval(t), derivative2.py.eval(t), derivative2.pz.eval(t)])

   # 3rd derivative
    derivative3 = derivative2.derivative()
    result.jerk = np.array([derivative3.px.eval(t), derivative3.py.eval(t), derivative3.pz.eval(t)])
    
    # 4th derivative
    derivative4 = derivative3.derivative()
    result.snap = np.array([derivative4.px.eval(t), derivative4.py.eval(t), derivative4.pz.eval(t)])

    thrust = result.acc + np.array([0, 0, 9.81]) # add gravity

    z_body = normalize(thrust)
    x_world = np.array([np.cos(result.yaw), np.sin(result.yaw), 0])
    y_body = normalize(np.cross(z_body, x_world))
    x_body = np.cross(y_body, z_body)

    jerk_orth_zbody = result.jerk - (np.dot(result.jerk, z_body) * z_body)
    h_w = jerk_orth_zbody / norm(thrust)

    result.omega = np.array([-np.dot(h_w, y_body), np.dot(h_w, x_body), z_body[2] * dyaw])

    # compute required roll/pitch angles
    result.pitch = np.arcsin(-x_body[2])
    result.roll = np.arctan2(y_body[2], z_body[2])

    return result


  # Matthias Faessler, Antonio Franchi, and Davide Scaramuzza:
  #     Differential Flatness of Quadrotor Dynamics Subject to Rotor Drag for Accurate Tracking of High-Speed Trajectories.

  def eval1(self, mass, t):
    result = TrajectoryOutput()
    # flat variables
    result.pos = np.array([self.px.eval(t), self.py.eval(t), self.pz.eval(t)])
    result.yaw = self.pyaw.eval(t)

    # 1st derivative
    derivative = self.derivative()
    result.vel = np.array([derivative.px.eval(t), derivative.py.eval(t), derivative.pz.eval(t)])
    dyaw = derivative.pyaw.eval(t)

    # 2nd derivative
    derivative2 = derivative.derivative()
    result.acc = np.array([derivative2.px.eval(t), derivative2.py.eval(t), derivative2.pz.eval(t)])

    # 3rd derivative
    derivative3 = derivative2.derivative()
    result.jerk = np.array([derivative3.px.eval(t), derivative3.py.eval(t), derivative3.pz.eval(t)])
    
    # 4th derivative
    derivative4 = derivative3.derivative()
    result.snap = np.array([derivative4.px.eval(t), derivative4.py.eval(t), derivative4.pz.eval(t)])
    
    x_c = np.array([np.cos(result.yaw), np.sin(result.yaw), 0]) 
    y_c = np.array([-np.sin(result.yaw), np.cos(result.yaw), 0]) 

    alpha = result.acc + np.array([0,0,9.81])
    beta  = alpha
    x_body = normalize(np.cross(y_c, alpha))
    y_body = normalize(np.cross(beta,  x_body))
    z_body = np.cross(x_body, y_body)
    # compute required roll/pitch angles
    result.pitch = np.arcsin(-x_body[2])
    result.roll = np.arctan2(y_body[2], z_body[2])

    c = np.dot(z_body, alpha)
    
    wx  = (-c * norm(np.cross(y_c, z_body)) * np.dot(y_body, result.jerk)) / (c * (c * norm(np.cross(y_c, z_body))))
    wy  = (np.dot(x_body, result.jerk))/(c)
    wz  = ((c * dyaw * np.dot(x_c, x_body)) + (np.dot(y_c, z_body) * np.dot(x_body, result.jerk))) / (c * norm(np.cross(y_c, z_body)))
    
    result.omega = np.array([wx, wy, wz])

    return result

  # A Comparative Study of Nonlinear MPC and Differential-Flatness-Based Control for Quadrotor Agile Flight
  # Sihao Sun, Angel Romero, Philipp Foehn, Elia Kaufmann and Davide Scaramuzza
  def eval2(self, mass, t):
    result = TrajectoryOutput()
    # flat variables
    result.pos = np.array([self.px.eval(t), self.py.eval(t), self.pz.eval(t)])
    result.yaw = self.pyaw.eval(t)

    # 1st derivative
    derivative = self.derivative()
    result.vel = np.array([derivative.px.eval(t), derivative.py.eval(t), derivative.pz.eval(t)])
    dyaw = derivative.pyaw.eval(t)

    # 2nd derivative
    derivative2 = derivative.derivative()
    result.acc = np.array([derivative2.px.eval(t), derivative2.py.eval(t), derivative2.pz.eval(t)])

    # 3rd derivative
    derivative3 = derivative2.derivative()
    result.jerk = np.array([derivative3.px.eval(t), derivative3.py.eval(t), derivative3.pz.eval(t)])
    
    # 4th derivative
    derivative4 = derivative3.derivative()
    result.snap = np.array([derivative4.px.eval(t), derivative4.py.eval(t), derivative4.pz.eval(t)])

    F = mass * (result.acc + np.array([0,0,9.81]))
   
    z_body = normalize(F)
    x_c = np.array([np.cos(result.yaw), np.sin(result.yaw), 0]) 
    y_body = normalize(np.cross(z_body, x_c))
    x_body = np.cross(y_body, z_body)
    z_world = np.array([0,0,1])
     
    # compute required roll/pitch angles
    result.pitch = np.arcsin(-x_body[2])
    result.roll = np.arctan2(y_body[2], z_body[2])

    thrust = norm(F)
    thrustDot = np.dot(mass*result.jerk, z_body)
    h_w = (mass/thrust)*(result.jerk - thrustDot*z_body)

    result.omega = np.array([-np.dot(h_w, y_body), np.dot(h_w, x_body), dyaw*np.dot(z_world, z_body)])
    
    return result

class Trajectory:
  def __init__(self):
    self.polynomials = None
    self.duration = None

  def loadcsv(self, filename):
    data = np.loadtxt(filename, delimiter=",", skiprows=1, usecols=range(33), ndmin=2)
    self.polynomials = [Polynomial4D(row[0], row[1:9], row[9:17], row[17:25], row[25:33]) for row in data]
    self.duration = np.sum(data[:,0])

  def savecsv(self, filename):
    data = np.empty((len(self.polynomials), 8*4+1))
    for i, p in enumerate(self.polynomials):
      data[i,0] = p.duration
      data[i,1:9] = p.px.p
      data[i,9:17] = p.py.p
      data[i,17:25] = p.pz.p
      data[i,25:33] = p.pyaw.p
    np.savetxt(filename, data, fmt="%.6f", delimiter=",", header="duration,x^0,x^1,x^2,x^3,x^4,x^5,x^6,x^7,y^0,y^1,y^2,y^3,y^4,y^5,y^6,y^7,z^0,z^1,z^2,z^3,z^4,z^5,z^6,z^7,yaw^0,yaw^1,yaw^2,yaw^3,yaw^4,yaw^5,yaw^6,yaw^7")
  
  def addmass(self, mass):
    self.mass = mass
  
  def stretchtime(self, factor):
    for p in self.polynomials:
      p.stretchtime(factor)
    self.duration *= factor

  def eval(self, t):
    assert t >= 0
    assert t <= self.duration

    current_t = 0.0
    for p in self.polynomials:
      if t < current_t + p.duration:
        return p.eval(t - current_t)
      current_t = current_t + p.duration

  # Matthias Faessler, Antonio Franchi, and Davide Scaramuzza:
  #     Differential Flatness of Quadrotor Dynamics Subject to Rotor Drag for Accurate Tracking of High-Speed Trajectories.
  def eval1(self, t):
    assert t >= 0
    assert t <= self.duration

    current_t = 0.0
    for p in self.polynomials:
      if t < current_t + p.duration:
        return p.eval1(self.mass, t - current_t)
      current_t = current_t + p.duration

  # A Comparative Study of Nonlinear MPC and Differential-Flatness-Based Control for Quadrotor Agile Flight
  # Sihao Sun, Angel Romero, Philipp Foehn, Elia Kaufmann and Davide Scaramuzza
  def eval2(self, t):
    assert t >= 0
    assert t <= self.duration

    current_t = 0.0
    for p in self.polynomials:
      if t < current_t + p.duration:
        return p.eval2(self.mass, t - current_t)
      current_t = current_t + p.duration

