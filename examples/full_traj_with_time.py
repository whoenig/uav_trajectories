import numpy as np



r = 1.0  # meter 
T = 10.0 # sec
omega = 10.0 # rad/sec


with open("timed_waypoints_helix_yaw.csv", "w") as f:
    f.write("t,x,y,z,yaw\n")

    for t in np.linspace(0, T, 100):
        f.write("{},{},{},{},{}\n".format(t, r * np.cos(t), r * np.sin(t), 0.2*t, omega*t))
