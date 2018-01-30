# trajectory
Estimate the trajectory of a LiDAR sensor based on multiple returns.

## Theory behind the code
When returns from a single pulse are detected, the sensor compute their positions as being in the center of the footprint and thus being all aligned. Because of that beheaviour, a line drawn between and beyond those returns must cross the sensor. Thus, several consecutive pulses emitted in a tight interval (e.g. 0.001 second) can be used to approximate an intersection point in the sky that would correspond to the sensor position given that the sensor carrier hasn't move much during this interval. A weighed least squares method using pseudoinverse gives a "close enough" approximation by minimising the squared sum of the distances between the intersection point and all the lines. Returns from a single pulse that are shortly spaced can be removed too prevent a large uncertainty on the intersection point due to precision error on coordinates of returns. In other words, two points spaced by 10 meters with 0.001 meter precision on XYZ coordinates will yield a more accurate line in the sky than two points spaced by 0.5 meter with the same precision on XYZ coordinates.

##### Written for
*Ministère des Forêts, de la Faune et des Parcs, Gouvernement du Québec*\
*Direction des inventaires forestiers*
