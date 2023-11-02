#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import argparse

import uav_trajectory

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("trajectory", type=str, help="CSV file containing trajectory")
  parser.add_argument("--dt", type=float, default=0.01, help="time steps in seconds")
  parser.add_argument("--mass", type=float, help="mass of the UAV")
  parser.add_argument("--stretchtime", type=float, help="stretch time factor (smaller means faster)")
  args = parser.parse_args()

  traj = uav_trajectory.Trajectory()
  traj.loadcsv(args.trajectory)
  pdfname = args.trajectory.replace('.csv','.pdf')
  f = PdfPages(pdfname)

  if args.mass:
    traj1 = uav_trajectory.Trajectory()
    traj1.loadcsv(args.trajectory)
    traj2 = uav_trajectory.Trajectory()
    traj2.loadcsv(args.trajectory)
    traj1.addmass(args.mass)
    traj2.addmass(args.mass)
    
  if args.stretchtime:
    traj.stretchtime(args.stretchtime)
    traj1.stretchtime(args.stretchtime)
    traj2.stretchtime(args.stretchtime)
    
  ts = np.arange(0, traj.duration, args.dt)
  evals = np.empty((len(ts), 21))
  for t, i in zip(ts, range(0, len(ts))):
    e = traj.eval(t)
    evals[i, 0:3]   = e.pos
    evals[i, 3:6]   = e.vel
    evals[i, 6:9]   = e.acc
    evals[i, 9:12]  = e.omega
    evals[i, 12]    = e.yaw
    evals[i, 13]    = e.roll
    evals[i, 14]    = e.pitch
    evals[i, 15:18] = e.jerk
    evals[i, 18::]  = e.snap
  
  full_traj = np.zeros((16, len(ts)))
  full_traj[0,:] = ts
  full_traj[1:4,:]  = np.transpose(evals[:,0:3])
  full_traj[4:7,:]  = np.transpose(evals[:,3:6])
  full_traj[7:10,:] = np.transpose(evals[:,6:9])
  full_traj[10:13,:] = np.transpose(evals[:,9:12])
  full_traj[13:16,:] = np.transpose(evals[:,12:15])

  filename = args.trajectory
  filename = filename.replace('.csv','_traj.csv')
  np.savetxt(filename, full_traj, delimiter=',')

  if args.mass:
    evals1 = np.empty((len(ts), 21))
    evals2 = np.empty((len(ts), 21))
    for t, i in zip(ts, range(0, len(ts))):
      e1 = traj1.eval1(t)
      e2 = traj2.eval2(t)
      evals1[i, 0:3]  = e1.pos
      evals1[i, 3:6]  = e1.vel
      evals1[i, 6:9]  = e1.acc
      evals1[i, 9:12] = e1.omega
      evals1[i, 12]   = e1.yaw
      evals1[i, 13]   = e1.roll
      evals1[i, 14]   = e1.pitch
      evals1[i, 15:18] = e1.jerk
      evals1[i, 18::]  = e1.snap
      evals2[i, 0:3]  = e2.pos
      evals2[i, 3:6]  = e2.vel
      evals2[i, 6:9]  = e2.acc
      evals2[i, 9:12] = e2.omega
      evals2[i, 12]   = e2.yaw
      evals2[i, 13]   = e2.roll
      evals2[i, 14]   = e2.pitch
      evals2[i, 15:18] = e2.jerk
      evals2[i, 18::]  = e2.snap
    velocity1 = np.linalg.norm(evals1[:,3:6], axis=1)
    velocity2 = np.linalg.norm(evals2[:,3:6], axis=1)
    acceleration1 = np.linalg.norm(evals1[:,6:9], axis=1)
    acceleration2 = np.linalg.norm(evals2[:,6:9], axis=1)
    omega1 = np.linalg.norm(evals1[:,9:12], axis=1)
    omega2 = np.linalg.norm(evals2[:,9:12], axis=1)
    
  velocity = np.linalg.norm(evals[:,3:6], axis=1)
  acceleration = np.linalg.norm(evals[:,6:9], axis=1)
  omega = np.linalg.norm(evals[:,9:12], axis=1)

  # print stats
  if args.mass:
    print("max speed (m/s): (W., Fr., Scr.)", np.max(velocity), np.max(velocity1), np.max(velocity2))
    print("max acceleration (m/s^2): (W., Fr., Scr.) ", np.max(acceleration), np.max(acceleration1), np.max(acceleration2))
    print("max omega (rad/s): (W., Fr., Scr.) ", np.max(omega), np.max(omega1), np.max(omega2))
    print("max roll  (deg): (W., Fr., Scr.) ", np.max(np.degrees(evals[:,13])), np.max(np.degrees(evals1[:,13])), np.max(np.degrees(evals2[:,13])))
    print("max pitch (deg): (W., Fr., Scr.) ", np.max(np.degrees(evals[:,14])), np.max(np.degrees(evals1[:,14])), np.max(np.degrees(evals2[:,14])))
    print("max yaw   (deg): (W., Fr., Scr.) ", np.max(np.degrees(evals[:,12])), np.max(np.degrees(evals1[:,12])), np.max(np.degrees(evals2[:,12])))
  else:
    print("max speed (m/s): ", np.max(velocity))
    print("max acceleration (m/s^2): ", np.max(acceleration))
    print("max omega (rad/s): ", np.max(omega))
    print("max roll (deg): ", np.max(np.degrees(evals[:,13])))
    print("max pitch (deg): ", np.max(np.degrees(evals[:,14])))

  # Create 3x1 sub plots
  gs = gridspec.GridSpec(6, 1)
  fig1 = plt.figure()

  ax1 = plt.subplot(gs[0:2, 0], projection='3d') # row 0
  ax1.plot(evals[:,0], evals[:,1], evals[:,2], lw=0.85)

  ax1 = plt.subplot(gs[2, 0]) # row 2
  ax1.plot(ts, velocity, lw=0.85)
  ax1.set_ylabel("vel")

  ax1 = plt.subplot(gs[3, 0]) # row 3
  ax1.plot(ts, acceleration, lw=0.85)
  ax1.set_ylabel("acc")

  ax1 = plt.subplot(gs[4, 0]) # row 4
  ax1.plot(ts, omega, lw=0.85)
  if args.mass:
    ax1.plot(ts, omega1, lw=0.85)
    ax1.plot(ts, omega2, lw=0.85)
  ax1.set_ylabel("w")

  ax1 = plt.subplot(gs[5, 0]) # row 5
  ax1.plot(ts, np.degrees(evals[:,12]), lw=0.85)
  ax1.set_ylabel("yaw ")

  if args.mass:
    gs = gridspec.GridSpec(3, 1)
    fig2 = plt.figure()
    plt.legend('omega_i')

    ax2 = plt.subplot(gs[0, 0])
    ax2.plot(ts, evals[:,9] ,lw=0.85)
    ax2.plot(ts, evals1[:,9] ,lw=0.85)
    ax2.plot(ts, evals2[:,9] ,lw=0.85)
    ax2.set_ylabel('omega_x')
    
    ax2 = plt.subplot(gs[1, 0])
    ax2.plot(ts, evals[:,10] , lw=0.85)
    ax2.plot(ts, evals1[:,10], lw=0.85)
    ax2.plot(ts, evals2[:,10], lw=0.85)
    ax2.set_ylabel('omega_y')

    ax2 = plt.subplot(gs[2, 0])
    ax2.plot(ts, evals[:,11] , label='W.', lw=0.85)
    ax2.plot(ts, evals1[:,11], label='Fr.',lw=0.85)
    ax2.plot(ts, evals2[:,11], label='Sc.',lw=0.85)
    ax2.set_ylabel('omega_z')
    ax2.legend(loc='lower center', ncol=3)

  fig1.savefig(f, format='pdf', bbox_inches='tight')
  fig2.savefig(f, format='pdf', bbox_inches='tight')
  f.close()
  plt.show()


