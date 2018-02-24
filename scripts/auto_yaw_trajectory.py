#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import argparse
# import scipy.interpolate
import scipy.optimize

import uav_trajectory

def func(coefficients, tss, yawss):
  result = 0
  for ts, yaws, i in zip(tss, yawss, range(0, len(tss))):
    yaws_output = np.polyval(coefficients[i*8:(i+1)*8], ts)
    result += np.sum((yaws - yaws_output) ** 2)
    # print(yaws_output)
  return result

def func_eq_constraint_val(coefficients, i, tss, yawss):
  result = 0
  end_val = np.polyval(coefficients[(i-1)*8:i*8], tss[i-1][-1])
  start_val = np.polyval(coefficients[i*8:(i+1)*8], tss[i][0])
  # print(i, end_val, start_val)
  return end_val - start_val

def func_eq_constraint_der(coefficients, i, tss, yawss):
  result = 0
  last_der = np.polyder(coefficients[(i-1)*8:i*8])
  this_der = np.polyder(coefficients[i*8:(i+1)*8])

  end_val = np.polyval(last_der, tss[i-1][-1])
  start_val = np.polyval(this_der, tss[i][0])
  return end_val - start_val

def func_eq_constraint_der_value(coefficients, i, t, desired_value):
  result = 0
  der = np.polyder(coefficients[i*8:(i+1)*8])

  value = np.polyval(der, t)
  return value - desired_value

# def func_eq_constraint(coefficients, tss, yawss):
#   result = 0
#   last_derivative = None
#   for ts, yaws, i in zip(tss, yawss, range(0, len(tss))):
#     derivative = np.polyder(coefficients[i*8:(i+1)*8])
#     if last_derivative is not None:
#       result += np.polyval(derivative, 0) - last_derivative
#     last_derivative = np.polyval(derivative, tss[-1])


  # # apply coefficients to trajectory
  # for i,p in enumerate(traj.polynomials):
  #   p.pyaw.p = coefficients[i*8:(i+1)*8]
  # # evaluate at each timestep and compute the sum of squared differences
  # result = 0
  # for t,yaw in zip(ts,yaws):
  #   e = traj.eval(t)
  #   result += (e.yaw - yaw) ** 2
  # return result

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("trajectory", type=str, help="CSV file containing trajectory")
  parser.add_argument("output", type=str, help="CSV file containing trajectory with updated yaw")
  parser.add_argument("--num", type=int, default=20, help="number of sampled points per trajectory segment")
  args = parser.parse_args()

  traj = uav_trajectory.Trajectory()
  traj.loadcsv(args.trajectory)

  tss = []
  yawss = []
  for p in traj.polynomials:
    ts = np.linspace(0, p.duration, args.num) #np.arange(0, p.duration, args.dt)
    evals = np.empty((len(ts), 15))
    for t, i in zip(ts, range(0, len(ts))):
      e = p.eval(t)
      evals[i, 0:3]  = e.pos
      evals[i, 3:6]  = e.vel
      evals[i, 6:9]  = e.acc
      evals[i, 9:12] = e.omega
      evals[i, 12]   = e.yaw
      evals[i, 13]   = e.roll
      evals[i, 14]   = e.pitch

    yaws = np.arctan2(evals[:,4], evals[:,3])
    tss.append(ts)
    yawss.append(yaws)

  x0 = np.zeros(len(traj.polynomials) * 8)
  print(x0)

  constraints = []
  for i in range(1, len(tss)):
    constraints.append({'type': 'eq', 'fun': func_eq_constraint_val, 'args': (i, tss, yawss)})
    constraints.append({'type': 'eq', 'fun': func_eq_constraint_der, 'args': (i, tss, yawss)})

  # zero derivative at the beginning and end
  constraints.append({'type': 'eq', 'fun': func_eq_constraint_der_value, 'args': (0, tss[0][0], 0)})
  constraints.append({'type': 'eq', 'fun': func_eq_constraint_der_value, 'args': (len(tss)-1, tss[-1][-1], 0)})


  res = scipy.optimize.minimize(func, x0, (tss, yawss), method="SLSQP", options={"maxiter": 1000}, 
    constraints=constraints
    )
  print(res)

  for i,p in enumerate(traj.polynomials):
    result = res.x[i*8:(i+1)*8]
    p.pyaw.p = np.array(result[::-1])

  traj.savecsv(args.output)

  # current_t = 0.0
  # for p in traj.polynomials:
  #   ts = np.arange(0, p.duration, args.dt)
  #   evals = np.empty((len(ts), 15))
  #   for t, i in zip(ts, range(0, len(ts))):
  #     e = p.eval(t)
  #     evals[i, 0:3]  = e.pos
  #     evals[i, 3:6]  = e.vel
  #     evals[i, 6:9]  = e.acc
  #     evals[i, 9:12] = e.omega
  #     evals[i, 12]   = e.yaw
  #     evals[i, 13]   = e.roll
  #     evals[i, 14]   = e.pitch

  #   # velocity = np.linalg.norm(evals[:,3:6], axis=1)
  #   yaw = np.arctan2(evals[:,4], evals[:,3])
  #   print(yaw)
  #   # yaw[0] = 0
  #   # yaw[-1] = 0
  #   result = np.polyfit(ts, yaw, deg=7)
  #   # p.pyaw.p = np.append(result[::-1], [0,0,0])
  #   p.pyaw.p = np.array(result[::-1])
  #   # print(result)
  #   # result = scipy.interpolate.UnivariateSpline(ts, yaw, k=5)
  #   # print(result(ts))

  #   # p.pyaw.p = np.append(result.get_coeffs(), [0, 0])
  #   # print(p.pyaw.p)
  #   # print(result.get_coeffs())
  # traj.savecsv(args.output)





  # print(yaw)



