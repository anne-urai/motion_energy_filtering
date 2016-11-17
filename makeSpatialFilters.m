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

function [kernel1, kernel2, cfg] = makeSpatialFilters(cfg)

% set defaults - from Kiani et al. 2008
% sigma_g, see email Klaus 29/09/2016
if ~isfield(cfg, 'sigma_c'),    cfg.sigma_c    = 0.35;                  end
if ~isfield(cfg, 'sigma_g'),    cfg.sigma_g    = 0.35;                  end
if ~isfield(cfg, 'x_range'),    cfg.x_range    = [-1 1];                end
if ~isfield(cfg, 'ppd'),        error('pixels per degree required!');   end

% formula from Kiani et al. 2008
alpha   = @(x) atand(x ./ cfg.sigma_c);
f1      = @(x,y) (cosd(alpha(x)).^4 .* cosd(4*alpha(x)) .* exp(-(y.^2)/(2*cfg.sigma_g^2)))';
f2      = @(x,y) (cosd(alpha(x)).^4 .* sind(4*alpha(x)) .* exp(-(y.^2)/(2*cfg.sigma_g^2)))';

% number of points for the kernel in x and y 
% (looks complicated but just ensures an odd number of points so that the center is at zero)
cfg.x_n     = 2*floor(floor((cfg.x_range(2) - cfg.x_range(1)) * cfg.ppd) / 2) + 1;  

% make spatial grid
space   = linspace(cfg.x_range(1), cfg.x_range(2),cfg. x_n);
[x, y]  = meshgrid(space, space);

% get the kernels
kernel1 = f1(x,y);
kernel2 = f2(x,y);

end

