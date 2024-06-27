function HOLD_PLOT_VEC3(vec, str)
x = [0 vec(1)];
y = [0 vec(2)];
z = [0 vec(3)];
plot3(x,y,z, 'w--', 'LineWidth',1, 'MarkerSize',10);
plot3(vec(1),vec(2),vec(3), 'r.','MarkerSize',15);
text(vec(1),vec(2),vec(3),string(str))

