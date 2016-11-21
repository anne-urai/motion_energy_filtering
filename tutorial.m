
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai and Klaus Wimmer, 2016
% anne.urai@gmail.com / kwimmer@clinic.ub.es
%
% ======================================================= %
% MOTION ENERGY FILTERING as described by Adelson & Bergen (1985)
% ======================================================= %
%
% see also these excellent tutorials, mainly for 1D motion
% George Mather
%   http://www.georgemather.com/Model.html
% Simoncelli, Glimcher, Chichilnisky
%   http://mplab.ucsd.edu/~marni/CSHL_Tutorials/MotionEnergy.tar
% David Heeger
%   http://www.cns.nyu.edu/~david/courses/perceptionGrad/syllabus2015.html
%   http://www.cns.nyu.edu/~david/courses/matlabTutorials/motionTutorial.zip

clear all; close all; clc;

% ======================================================= %
% GET COORDINATES AND TURN INTO DOT MOVIE
% ======================================================= %

thetas = 0:45:360; me = nan(size(thetas));

for theta = 1:length(thetas),
    
    % simulate some dots
    display.frameRate = 60; % Hz
    display.width     = 41; % cm
    display.dist      = 50; % cm
    display.res.width = 1280; % pixels
    display.center    = [640 512];
    display.ppd       = deg2pix(display, 1);
    
    setup.nframes     = 30; %1 * display.frameRate; % shorter for testing quickly
    dots              = setupDots(display);
    dots.direction    = 0; % 0 = rightwards
    dots.coherence    = 1; % very strong coherence
    dots.nvar         = 1; % so speed we input is really the effective speed
    dots.lifetime     = setup.nframes;
    coord             = dots_limitedlifetime(setup, display, dots);
    
    % create a stimulus 3d movie from coordinates
    % coordinates must be nframes * x/y * nDots
    % motion will be turned into the downwards direction, easiest for filtering
    stim              = coord2stim(display, coord, -dots.direction + thetas(theta));
    
    % show movie - make sure the dots move to the right!
    % save2gif(stim);
    
    % show the spatiotemporal fourier representation along x axis
    % spectogram(squeeze(sum(stim,2))', display);
    % title(sprintf('%2.2f deg/s', spd / dots.nvar * cfg.frameRate / cfg.ppd));
    
    % otherwise, load your stimulus and display settings here
    
    % ======================================================= %
    % CREATE SPATIAL AND TEMPORAL FILTERS
    % ======================================================= %
    
    % temporal range of the filter
    cfg.frameRate  = display.frameRate;
    cfg.ppd        = display.ppd;
    
    % k = 60, from Kiani et al. 2008
    cfg.k = 60;
    
    % adjust spatial filters to match the speed in the dots
    effectiveSpeed = pix2deg(display, dots.speed) ./ dots.nvar;
    
    % if we only have coordinates, can try to reconstruct the generative speed
    % warning: if multiple sets of frames were interleaved, this doesn't work!
    % computedSpeed = computeEffectiveSpeed(coord, display);
    % assert(abs(effectiveSpeed - computedSpeed) < 0.1, 'could not compute generative speed!');
    
    % Kiani et al. 2008 has a speed of 2.84 deg/s and used sigma_c and sigma_g
    % as 0.35 (not explicitly mentioned in their paper). To optimally scale the
    % filters for the speed in our dots, multiply the spatial parameters
    cfg.sigma_c = 0.35 * (effectiveSpeed / 2.84);
    cfg.sigma_g = 0.35 * (effectiveSpeed / 2.84);
    
    % equations exactly as in Kiani et al. 2008
    [f1, f2] = makeSpatialFilters(cfg);
    [g1, g2] = makeTemporalFilters(cfg);
    
    % ======================================================= %
    % FILTER TO OBTAIN MOTION ENERGY
    % ======================================================= %
    
    motionenergy = applyFilters(stim, f1, f2, g1, g2);
    me(theta) = sum(sum(sum(motionenergy))); % sum over space and time
end

polar(deg2rad(thetas), me);
polar(deg2rad(thetas), ones(size(thetas)));
