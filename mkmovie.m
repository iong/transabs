g = h5info('t.h5');
N = size(g.Groups);
F(N) = struct('cdata',[],'colormap',[]);
boxmax=30;
v = VideoWriter('test.avi');
open(v);
for i=1:N
    q = h5read('t.h5', sprintf('/%d/q', i));
    x = h5read('t.h5', sprintf('/%d/x', i));
    y = h5read('t.h5', sprintf('/%d/y', i));
    z = h5read('t.h5', sprintf('/%d/z', i));
    
    scatter3(x, y, z, 10*abs(q), q);
    axis([-boxmax boxmax -boxmax boxmax -boxmax boxmax]);
    F(i) = getframe;
    writeVideo(v, F(i));
end
close(v);
    