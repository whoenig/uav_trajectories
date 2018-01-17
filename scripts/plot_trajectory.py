#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import argparse

import uav_trajectory

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("trajectory", type=str, help="CSV file containing trajectory")
  parser.add_argument("mass", type=float, help="mass of UAV [kg]")
  args = parser.parse_args()

  traj = uav_trajectory.Trajectory()
  traj.loadcsv(args.trajectory)

  ts = np.arange(0, traj.duration, 0.01)
  evals = np.empty((len(ts), 13))
  for t, i in zip(ts, range(0, len(ts))):
    e = traj.eval(args.mass, t)
    evals[i, 0:3]  = e.pos
    evals[i, 3:6]  = e.vel
    evals[i, 6:9]  = e.acc
    evals[i, 9:12] = e.omega
    evals[i, 12]   = e.yaw

  # Create 3x1 sub plots
  gs = gridspec.GridSpec(6, 1)
  fig = plt.figure()

  ax = plt.subplot(gs[0:2, 0], projection='3d') # row 0
  ax.plot(evals[:,0], evals[:,1], evals[:,2])

  ax = plt.subplot(gs[2, 0]) # row 2
  ax.plot(ts, np.linalg.norm(evals[:,3:6], axis=1))
  ax.set_ylabel("velocity [m/s]")

  ax = plt.subplot(gs[3, 0]) # row 3
  ax.plot(ts, np.linalg.norm(evals[:,6:9], axis=1))
  ax.set_ylabel("acceleration [m/s^2]")

  ax = plt.subplot(gs[4, 0]) # row 4
  ax.plot(ts, np.linalg.norm(evals[:,9:12], axis=1))
  ax.set_ylabel("omega [rad/s]")

  ax = plt.subplot(gs[5, 0]) # row 5
  ax.plot(ts, evals[:,12])
  ax.set_ylabel("yaw [rad]")

  plt.show()
