from ubuntu:16.04

RUN apt-get update && apt-get install -y git build-essential libnlopt-dev libgoogle-glog-dev cmake libeigen3-dev libboost-.*-dev

ADD ./ /uav_trajectories/

RUN mkdir /uav_trajectories/build
WORKDIR /uav_trajectories/build

RUN cmake ..
RUN make

ENTRYPOINT [ "/uav_trajectories/build/genTrajectory" ]
