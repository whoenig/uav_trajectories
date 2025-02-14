#!/usr/bin/env python

import numpy as np
import rowan
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import argparse

import uav_trajectory

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("trajectory", type=str, help="CSV file containing trajectory")
  parser.add_argument("--stretchtime", type=float, help="stretch time factor (smaller means faster)")
  args = parser.parse_args()

  # Crazyflie 2.1 params
  mass = 0.034
  inertia = np.diag([16.571710e-6, 16.655602e-6, 29.261652e-6])
  arm_length = 0.046  # m
  arm = 0.707106781 * arm_length
  t2t = 0.006  # thrust-to-torque ratio
  allocation_matrix = np.array([
      [1, 1, 1, 1],
      [-arm, -arm, arm, arm],
      [-arm, arm, arm, -arm],
      [-t2t, t2t, -t2t, t2t]
      ])

  traj = uav_trajectory.Trajectory()
  traj.loadcsv(args.trajectory)

  if args.stretchtime:
    traj.stretchtime(args.stretchtime)

  ts = np.arange(0, traj.duration, 0.01)
  evals = np.empty((len(ts), 24))

  for t, i in zip(ts, range(0, len(ts))):
    e = traj.eval(t, mass, inertia)
    evals[i, 0:3]  = e.pos
    evals[i, 3:6]  = e.vel
    evals[i, 6:10] = rowan.from_matrix(e.rotation)
    evals[i, 10:13] = e.omega

    evals[i, 13:14] = e.force
    evals[i, 14:17] = e.torque

    evals[i, 17:21] = np.linalg.inv(allocation_matrix) @ evals[i,13:17]

    evals[i, 21:24] = e.acc

  velocity = np.linalg.norm(evals[:,3:6], axis=1)
  acceleration = np.linalg.norm(evals[:,21:24], axis=1)
  omega = np.linalg.norm(evals[:,10:13], axis=1)
  rpy = rowan.to_euler(rowan.normalize(evals[:,6:10]), "xyz", axis_type="extrinsic")


  # print stats
  print("max speed (m/s): ", np.max(velocity))
  print("max acceleration (m/s^2): ", np.max(acceleration))
  print("max omega (rad/s): ", np.max(omega))
  print("max roll (deg): ", np.max(np.degrees(rpy[:,0])))
  print("max pitch (deg): ", np.max(np.degrees(rpy[:,1])))

  # Create 3x1 sub plots
  gs = gridspec.GridSpec(6, 1)
  fig = plt.figure()

  ax = plt.subplot(gs[0:2, 0], projection='3d') # row 0
  ax.plot(evals[:,0], evals[:,1], evals[:,2])

  ax = plt.subplot(gs[2, 0]) # row 2
  ax.plot(ts, velocity)
  ax.set_ylabel("velocity [m/s]")

  ax = plt.subplot(gs[3, 0]) # row 3
  ax.plot(ts, acceleration)
  ax.set_ylabel("acceleration [m/s^2]")

  ax = plt.subplot(gs[4, 0]) # row 4
  ax.plot(ts, omega)
  ax.set_ylabel("omega [rad/s]")

  ax = plt.subplot(gs[5, 0]) # row 5
  ax.plot(ts, np.degrees(rpy[:,2]))
  ax.set_ylabel("yaw [deg]")

  # ax = plt.subplot(gs[6, 0]) # row 5
  # ax.plot(ts, np.degrees(evals[:,13]))
  # ax.set_ylabel("roll [deg]")

  # ax = plt.subplot(gs[7, 0]) # row 5
  # ax.plot(ts, np.degrees(evals[:,14]))
  # ax.set_ylabel("pitch [deg]")

  # plot the states
  fig, ax = plt.subplots(4, 3, sharex='all', squeeze=False)

  for k, axis in enumerate(["x", "y", "z"]):
      ax[0,k].plot(ts, evals[:,k])
      ax[0,0].set_ylabel("position [m]")

      ax[1,k].plot(ts, evals[:,3+k])
      ax[1,0].set_ylabel("velocity [m/s]")

      ax[2,k].plot(ts, np.degrees(rpy[:,k]))
      ax[2,0].set_ylabel("rotation [deg]")

      ax[3,k].plot(ts, np.degrees(evals[:,10+k]))
      ax[3,0].set_ylabel("angular velocity [deg/s]")

  # plot the force/torque
  fig, ax = plt.subplots(2, 2, sharex='all', squeeze=False)
  ax[0,1].plot(ts, evals[:,17])
  ax[0,1].set_ylabel("M1 [N]")
  ax[1,1].plot(ts, evals[:,18])
  ax[1,1].set_ylabel("M2 [N]")
  ax[1,0].plot(ts, evals[:,19])
  ax[1,0].set_ylabel("M3 [N]")
  ax[0,0].plot(ts, evals[:,20])
  ax[0,0].set_ylabel("M4 [N]")

  # print stats
  print("max motor force (N): ", np.max(evals[:,17:21]))

  plt.show()