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

function speed = computeEffectiveSpeed(coord, display)
% computes temporal tuning parameter k, based on either the speed of the
% dots (in degrees/s) or a set of dot coordinates 

% compute the displacement per dot between each pair of successive frames
moved = nan(size(coord, 3), size(coord, 1) -1);
for d = 1:size(coord, 3),
    for f = 1:size(coord, 1)-1,
        % displacement in x and y direction
        vecX = diff(coord([f f+1], 1, d));
        vecY = diff(coord([f f+1], 2, d));
        % compute the resulting vector
        vecResult = sqrt(vecX.^2 + vecY.^2);
        moved(d, f) = vecResult;
    end
end
moved = moved(:);

% now compute the speed from this displacement
speed        = pix2deg(display, moved(:)) .* display.frameRate; 
speed        = median(speed); % most prevalent speed in this display

end