function plot_phase_portait_pre_and_post_training(REC, imin, icrit)
    pre_time = 1:imin;
    post_time = 1:icrit;
    %% pre training
    for j = 1:1:5
        figure;
        plot(REC(pre_time,j), REC(pre_time,j+5))
        xlabel('Membrane Potential')
        ylabel('Adaption Variable') 
        title(strcat('Pre Learning Phase Space ', num2str(j)))
    end
    
    %% post training
    for j = 1:1:5
        figure;
        plot(REC(post_time,j), REC(post_time,j+5)) 
        xlabel('Membrane Potential')
        ylabel('Adaption Variable') 
        title(strcat('Post Learning Phase Space ', num2str(j)))
    end
    
end