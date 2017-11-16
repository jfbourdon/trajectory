# trajectory
Estimate the trajectory of a LiDAR sensor based on multiple returns

##### Written by
*Ministère des Forêts, de la Faune et des Parcs, Gouvernement du Québec*\
*Direction des inventaires forestiers*

## Theory behind the code
As all returns from a single pulse all aligned, the direction vector formed by thoses returns should cross the position of the sensor. Therefore, by taking several consecutive pulses, it is possible to evaluate the intersection in 3D space of all those direction vectors and thus estimate the trajectory of the LiDAR sensor for a given time interval.
