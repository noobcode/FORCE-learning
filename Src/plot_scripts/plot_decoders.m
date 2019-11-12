function plot_decoders(RECB, dt, nt)
    figure;
    plot((1:1:nt)*dt/1000, RECB)
    xlabel('Time (s)')
    ylabel('Decoders $\phi(t)$', 'Interpreter', 'LaTeX')
    title('Decoders')

end