import numpy as np

numRobots = 5

r = 0.3
height = 0.7
w = 2 * np.pi / numRobots
T = 2 * 2 * np.pi / w

# horizontal circles
for i in range(0, numRobots):
	phase = 2 * np.pi / numRobots * i

	with open("timed_waypoints_circle{}.csv".format(i), "w") as f:
		f.write("t,x,y,z\n")

		for t in np.linspace(0, T, 100):
			f.write("{},{},{},{}\n".format(t, r * np.cos(w * t + phase), r * np.sin(w * t + phase), height))
