# motionEnergy

Implements motion energy filtering as described in Adelson & Bergen (1985), on 3D moving dot stimuli as done in Kiani et al. (2008) and others.

See also these excellent tutorials on the 2D equivalent of these filters
* http://www.georgemather.com/Model.html (George Mather)
* http://mplab.ucsd.edu/~marni/CSHL_Tutorials/MotionEnergy.tar (Simoncelli, Glimcher, Chichilnisky)

INPUT
- stimulus, a 3d matrix containing the x-y-t stimulus movie
- display, structure containing the following fields:
                      - display.frameRate (in Hz)
                      - display.ppd (pixels/degree) OR all of the following:
                      - display.width (in cm of physical screen)
                      - display.res.width (in pixels)
                      - display.dist (between eyes and screen, in cm)
- theta, the angles at which motion energy should be computed
  (rightwards = 0, counterclockwise).
- plotme, generate graphical output (default = false);

OUTPUT
- motionenergy, length(theta) x size(stimulus, 3) motionenergy filter output
filters (optional)

Anne Urai, 2016
anne.urai@gmail.com / anneurai.net
