# Stochastic Background Noise Analysis
Laser Interferometer Space Antenna (LISA) Data Processing in Julia

## Description
Using a 3-week simulated sample from the planned LISA mission, NoiseData.dat, this script analyzses the noise from the detector in the data. Using a stochastic background, it searches for the noise levels from the gravitational wave background.

The stochastic background can be modelled by S<sub>gw</sub> = (3H<sup>2</sup><sub>0</sub>/4*pi<sup>2</sup>f<sup>3</sup>)*Ω<sub>gw</sub>, where H<sub>0</sub> is the hubble constant and Ω<sub>gw</sub> is the gravitational wave background level.

The expression of noise levels from the gravitational wave background combined with the detector response is XX<sub>gw</sub> = 4 * sin<sup>2</sup>(f/f<sub>star</sub>) * LowXX * S<sub>gw</sub>.

The script searches over both the gravitational background level and the acceleration noise and position noise from the detector.
