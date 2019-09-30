X_start = [0.51];
Y_start = [0.09];
Z_start = [0.71];

X_1 = [0.35];
Y_1 = [-0.2];
Z_1 = [0.1];

X_2 = [0.75];
Y_2 = [-0.2];
Z_2 = [0.1];

X_3 = [0.65];
Y_3 = [0.33];
Z_3 = [0.1];

X_4 = [0.35];
Y_4 = [0.33];
Z_4 = [0.1];

X_proj = [0.51];
Y_proj = [0.09];
Z_proj = [0.1];


figure;
scatter3(X_start,Y_start,Z_start, 50, 'filled')
hold on
scatter3(X_1,Y_1,Z_1, 50, 'LineWidth', 2)
scatter3(X_2,Y_2,Z_2, 50, 'LineWidth', 2)
scatter3(X_3,Y_3,Z_3, 50, 'LineWidth', 2)
scatter3(X_4,Y_4,Z_4, 50, 'LineWidth', 2)
scatter3(X_proj,Y_proj,Z_proj, 50, 'x')
plot3([X_1,X_2], [Y_1, Y_2], [Z_1, Z_2],'--', 'color', 'k')
plot3([X_1,X_4], [Y_1, Y_4], [Z_1, Z_4],'--', 'color', 'k')
plot3([X_2,X_3], [Y_2, Y_3], [Z_2, Z_3],'--', 'color', 'k')
plot3([X_3,X_4], [Y_3, Y_4], [Z_3, Z_4],'--', 'color', 'k')
plot3([X_start,X_proj], [Y_start, Y_proj], [Z_start, Z_proj],'--', 'color', 'b')
fill3([X_1, X_2, X_3, X_4], [Y_1, Y_2, Y_3, Y_4] , [Z_1, Z_2, Z_3, Z_4], 'green', 'FaceAlpha', 0.5)

legend('Start Point', 'End Point 1', 'End Point 2', 'End Point 3', 'End Point 4')
title('End-effector start and target positions')

zlabel('z')
xlabel('x')
ylabel('y')

