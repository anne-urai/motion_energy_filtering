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

function stim = coord2stim(display, coord, theta)
% transform coordinates into a 3d stimulus movie
% theta is the angle of dot motion, which can be used to align all
% motion in the zero direction (easier for filtering)

nFrames   = size(coord, 1);
nDots     = size(coord, 3);
stim      = zeros(display.center(1)*2, display.center(2)*2, nFrames);

% how big is the circle?
dotsRadius = max(abs(ceil(coord(:))));

for f = 1:nFrames,
    % rotate the whole stimulus movie towards the zero direction
    [th,r] = cart2pol(squeeze(coord(f, 1, :)), ...
        squeeze(coord(f, 2, :)));
    th = th + deg2rad(theta); % rotate by the indicated amount, theta in radians
    [xpos, ypos] = pol2cart(th, r); % convert back to cartesian coords
    
    % move towards the center, can't have negative screen indices
    xpos = round(xpos + display.center(1));
    ypos = round(ypos + display.center(2));
    
    % put those coordinates in the stim representation matrix
    for i = 1:nDots, % for each dot
        stim(xpos(i), ypos(i), f) = 1; % put a pixel in the matrix
    end
end

% first, remove non-informative border and make image square
x2use = display.center(1)-dotsRadius : display.center(1)+dotsRadius;
y2use = display.center(2)-dotsRadius : display.center(2)+dotsRadius;
stim  = stim(x2use, y2use, :);

% zero-padding for convolution, avoid edge artefacts in spatial dimension
% full convolution = [ma + mb - 1]
% assert 'same' == 'valid', see applyfilters c
stimpad    = roundn(2*floor(floor((2) * display.ppd) / 2) + 1, 2);
stim       = padarray(stim, [stimpad stimpad], 0, 'both');

end

% % rotate x and y coordinates towards zero direction
% % this does the same as the code above, but less clear to read
% xpos = ceil( squeeze(coord(f, 1, :)) * cosd(theta) ...
%     - squeeze(coord(f, 2, :)) * sind(theta) ...
%     + display.center(1) );
% ypos = ceil( squeeze(coord(f, 1, :)) * sind(theta) ...
%     + squeeze(coord(f, 2, :)) * cosd(theta) ...
%     + display.center(2));
