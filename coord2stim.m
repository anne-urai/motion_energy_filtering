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

nFrames   = size(coord, 1);
nDots     = size(coord, 3);
stim      = zeros(display.center(1)*2, display.center(2)*2, nFrames);

% how big is the circle?
dotsRadius = max(abs(ceil(coord(:))));
stimpad    = 50; % padding for convolution

% size of the screen in pixels
stimsize  = 2*dotsRadius+1+stimpad;

for f = 1:nFrames,
    % rotate x and y coordinates
    xpos = ceil( squeeze(coord(f, 1, :)) * cosd(theta) ...
        - squeeze(coord(f, 2, :)) * sind(theta) ...
        + display.center(1) );
    ypos = ceil( squeeze(coord(f, 1, :)) * sind(theta) ...
        + squeeze(coord(f, 2, :)) * cosd(theta) ...
        + display.center(2));
    
    % put those coordinates in the stim representation matrix
    for i = 1:nDots, % for each dot
        stim(xpos(i), ypos(i), f) = 1; % put a pixel in the matrix
    end
end

% keep smallest bit of the screen, minus padding (avoid filter artefacts)
x2use = display.center(1)-dotsRadius-stimpad/2 : display.center(1)+dotsRadius+stimpad/2;
y2use = display.center(2)-dotsRadius-stimpad/2 : display.center(2)+dotsRadius+stimpad/2;
stim  = stim(x2use, y2use, :);
assert(size(stim, 1) == stimsize, 'size mismatch');
assert(size(stim, 2) == stimsize, 'size mismatch');


end
