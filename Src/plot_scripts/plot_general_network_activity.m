function plot_general_network_activity(y, dt, nt)
    figure;
    plot(dt*(1:1:nt)/1000, y(2:end), 'LineWidth', 2, 'Color', 'black')
    xlabel('Time (s)')
    ylabel('$y(t)$','Interpreter','LaTeX')
    title('General Network Activity')
    xlim([0,nt*dt/1000]);
end