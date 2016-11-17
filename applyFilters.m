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

function energy = applyFilters(stim, f1, f2, g1, g2)
% apply filtering according to Adelson & Bergen figure 18b
% composite filters!
tic;
stimsize = size(stim, 1);
nFrames  = size(stim, 3);

% spatial filtering
A0 = zeros(stimsize, stimsize, nFrames);
B0 = zeros(stimsize, stimsize, nFrames);
for f = 1:nFrames,
    A0(:,:,f) = conv2(stim(:,:,f), f1, 'same');
    B0(:,:,f) = conv2(stim(:,:,f), f2, 'same');
end

% temporal filtering
A0_     = reshape(A0, stimsize * stimsize, []);
B0_     = reshape(B0, stimsize * stimsize, []);

A1_     = filter(g1, 1, A0_, [], 2);
A1      = reshape(A1_, stimsize, stimsize, []);

A2_     = filter(g2, 1, A0_, [], 2);
A2      = reshape(A2_, stimsize, stimsize, []);

B1_     = filter(g1, 1, B0_, [], 2);
B1      = reshape(B1_, stimsize, stimsize, []);

B2_     = filter(g2, 1, B0_, [], 2);
B2      = reshape(B2_, stimsize, stimsize, []);

% Adelson and Bergen, fig 18b
% opponentenergy    = 4*(A1 .* B2 - A2 .* B1); % this is the opponent energy
% energy_left       = (A1-B2).^2 + (A2+B1).^2;
% energy_right      = (A1+B2).^2 + (A2-B1).^2;
% energy_opponent   = energy_right - energy_left; % opponent energy

energy              = A1 .* B2 - A2 .* B1;

% energy_static     = A1.^2 + A2.^2 + B1.^2 + B2.^2;
% velocity          = energy_opponent ./ energy_static;
% assert(~any(isnan(velocity)));
toc;

end