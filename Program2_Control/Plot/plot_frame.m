function plot_frame(p,R,l,w)

% Vectors
vecx = [p p+l*R(:,1)];
vecy = [p p+l*R(:,2)];
vecz = [p p+l*R(:,3)];

% Plot
plot3( [vecx(1,1) vecx(1,2)], [vecx(2,1) vecx(2,2)], [vecx(3,1) vecx(3,2)], 'r', 'LineWidth', w); hold on;
plot3( [vecy(1,1) vecy(1,2)], [vecy(2,1) vecy(2,2)], [vecy(3,1) vecy(3,2)], 'g', 'LineWidth', w); hold on;
plot3( [vecz(1,1) vecz(1,2)], [vecz(2,1) vecz(2,2)], [vecz(3,1) vecz(3,2)], 'b', 'LineWidth', w); hold on;

end

