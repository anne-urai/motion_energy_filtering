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

function spectogram(inp, cfg)
% FFT spectrum of stimulus
% http://www.gaussianwaves.com/2014/07/how-to-plot-fft-using-matlab-fft-of-basic-signals-sine-and-cosine-waves/

% get the power spectrum
X       = fftshift(abs(fftn(inp)));

% compute correct axes
t   = (0:size(inp,2)) * 1/cfg.frameRate;
x   = ((1:size(inp,1)) - round(size(inp, 1) / 2)) * 1 / cfg.ppd;
f   = 1./(t(2)-t(1))*((1:numel(t))-(numel(t)/2+1))/numel(t);
x   = 1./(x(2)-x(1))*((1:numel(x))-(numel(x)/2+1))/numel(x);

colormap(1-gray);
imagesc(x, f, squeeze(X));
xlabel('Spatial frequency (1/deg)');
ylabel('Temporal frequency (1/s)');
axis([-5 5 -10 10]);
axis ij; % axis image additionally makes square;

end
