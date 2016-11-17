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

function [kernel1, kernel2, cfg] = makeTemporalFilters(cfg)

% set defaults
if ~isfield(cfg, 'k'),          cfg.k          = 60;              end
if ~isfield(cfg, 't_range'),    cfg.t_range    = [0 0.3];         end
if ~isfield(cfg, 'nslow'),      cfg.nslow      = 3;               end
if ~isfield(cfg, 'nfast'),      cfg.nfast      = 5;               end
if ~isfield(cfg, 'frameRate'),  error('frameRate required!');     end

% formula's from Kiani et al. 2008
g1 = @(t) (cfg.k*t).^cfg.nslow .* exp(-cfg.k*t) .* (1/factorial(cfg.nslow) - ((cfg.k*t).^2) / factorial(cfg.nslow + 2));
g2 = @(t) (cfg.k*t).^cfg.nfast .* exp(-cfg.k*t) .* (1/factorial(cfg.nfast) - ((cfg.k*t).^2) / factorial(cfg.nfast + 2));

% number of points for the kernel in t 
% (looks complicated but just ensures an odd number of points so that the center is at zero)
t_n = ceil((cfg.t_range(2) - cfg.t_range(1)) * cfg.frameRate);  
t   = linspace(cfg.t_range(1), cfg.t_range(2), t_n)';

% get the kernels
kernel1 = g1(t);
kernel2 = g2(t);

end

