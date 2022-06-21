import numpy as np

omega = 10 #rad/sec
T = 10 # sec

with open("timed_waypoints_yaw.csv", "w") as f:
    f.write("t,x,y,z,yaw\n")

    for t in np.linspace(0, T, 100):
        f.write("{},{},{},{},{}\n".format(t, 0, 0, 0, omega*t))
