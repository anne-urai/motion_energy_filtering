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

function save2gif(x)
% from https://sukhbinder.wordpress.com/2012/10/09/gif-animation-in-matlab-in-3-steps/
% makes a random dot gif

delete('imagefile.gif');
colormap(flipud(gray));
% loop over frames
for k = 1:size(x,3),
    contourf(x(:,:,k));
    axis square;
    f = getframe;
    [im, cm] = rgb2ind(f.cdata, 2, 'nodither');
    
    % save
    if k == 1;
        imwrite(im, cm, 'imagefile.gif', 'gif','DelayTime', 1/60);
    else
        imwrite(im, cm, 'imagefile.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 1/60);
    end
end

end
