function plot_raster(tspike, N, ns, dt, nt, imin, icrit)
    figure;
    min_ = 0;
    max_ = 100;
    hold on
    patch([imin icrit icrit imin]*dt/1000, [min_ min_ max_ max_], [0.5,0.5,0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
    plot(tspike(1:ns,2)/1000, tspike(1:ns,1),'k.')
    hold off
    xlabel('Time (s)')
    ylabel('Neuron Index')
    title('Raster plot')
    ylim([0,100])
    xlim([0,nt*dt/1000]);
end