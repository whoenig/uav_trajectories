import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import argparse
import os
from matplotlib.backends.backend_pdf import PdfPages
import uav_trajectory

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--traj", type=str, help="CSV file containing polynomials")
    parser.add_argument("--output", type=str, help="CSV file containing pos, vel, acc, jerk, snap")
    parser.add_argument("--dt", default=0.01, type=float, help="CSV file containing polynomials")
    parser.add_argument("--stretchtime", default=1.0, type=float, help="stretch time factor (smaller means faster)")
    args = parser.parse_args()

    traj = uav_trajectory.Trajectory()
    traj.loadcsv(args.traj)


    traj.stretchtime(args.stretchtime)


    ts = np.arange(0, traj.duration, args.dt)
    # t, pos, vel, acc, jerk, snap: 1 + 15 
    evals = np.empty((len(ts), 15 + 1))  # Additional column for 't'
    with open(args.output, "w") as f:
        for t, i in zip(ts, range(0, len(ts))):
            e = traj.eval(t)
            evals[i, 0] = t
            # pos[2] = 0
            evals[i, 1:4]  = e.pos
            evals[i, 4:7]  = e.vel
            evals[i, 7:10] = e.acc
            evals[i, 10:13] = e.jerk
            evals[i, 13:16] = e.snap
        header = "t,posx,posy,posz,velx,vely,velz,accx,accy,accz,jerkx,jerky,jerkz,snapx,snapy,snapz"

        evals_with_header = np.vstack((header.split(','), evals))
        np.savetxt(f, evals_with_header, delimiter=",", fmt='%s')



    # Extract position, velocity, and acceleration, jerk, snap data from evals
    pos = evals[:, 1:4]
    vel = evals[:, 4:7]
    acc = evals[:, 7:10]
    jerk = evals[:, 10:13]
    snap = evals[:, 13:16]

    # Time steps (using dt from arguments)
    time_steps = ts

    # Create a PDF file to save the plots
    pdf_filename = os.path.splitext(args.output)[0] + '.pdf'
    with PdfPages(pdf_filename) as pdf:  
        axes = ["x", "y", "z"]
        # Page 1: pos
        fig, ax = plt.subplots(3, 1, figsize=(8, 12))
        for i in range(3):
            ax[i].plot(time_steps, pos[:, i], label=f'pos[{i}]', color='b')
            ax[i].set_xlabel('Time (s)')
            ax[i].set_ylabel(f' {axes[i]}')
            ax[i].legend()
            ax[i].grid(True)
        fig.suptitle('Position')
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: vel
        fig, ax = plt.subplots(3, 1, figsize=(8, 12))
        for i in range(3):
            ax[i].plot(time_steps, vel[:, i], label=f'vel[{i}]', color='b')
            ax[i].set_xlabel('Time (s)')
            ax[i].set_ylabel(f' {axes[i]}')
            ax[i].legend()
            ax[i].grid(True)
        fig.suptitle('Velocity')
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3: acc
        fig, ax = plt.subplots(3, 1, figsize=(8, 12))
        for i in range(3):
            ax[i].plot(time_steps, acc[:, i], label=f'acc[{i}]', color='b')
            ax[i].set_xlabel('Time (s)')
            ax[i].set_ylabel(f' {axes[i]}')
            ax[i].legend()
            ax[i].grid(True)
        fig.suptitle('Acceleration')
        pdf.savefig(fig)
        plt.close(fig)


        # Page 4: jerk
        fig, ax = plt.subplots(3, 1, figsize=(8, 12))
        for i in range(3):
            ax[i].plot(time_steps, jerk[:, i], label=f'jerk[{i}]', color='b')
            ax[i].set_xlabel('Time (s)')
            ax[i].set_ylabel(f' {axes[i]}')
            ax[i].legend()
            ax[i].grid(True)
        fig.suptitle('Jerk')
        pdf.savefig(fig)
        plt.close(fig)

        # Page 5: snap
        fig, ax = plt.subplots(3, 1, figsize=(8, 12))
        for i in range(3):
            ax[i].plot(time_steps, snap[:, i], label=f'snap[{i}]', color='b')
            ax[i].set_xlabel('Time (s)')
            ax[i].set_ylabel(f' {axes[i]}')
            ax[i].legend()
            ax[i].grid(True)
        fig.suptitle('Snap')
        pdf.savefig(fig)
        plt.close(fig)

