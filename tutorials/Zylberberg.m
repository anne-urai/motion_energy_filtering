spacing = -1:.01:1;
[x, y] = meshgrid(spacing);

alpha = 0.56;
sx = 0.02;
sy = 0.97;
theta = [0 100 150]
for thistheta = theta;
    
    % create filters
    even = cos(2*pi*alpha*yprime) .* exp(-yprime.^2 ./ sy^2) .* exp(-xprime.^2 / sx^2);
    
    s(find(thistheta==theta)) = subplot(3,3, find(thistheta==theta));
    imagesc(even);
    xlabel('x'); ylabel('y'); title(thistheta)
end
