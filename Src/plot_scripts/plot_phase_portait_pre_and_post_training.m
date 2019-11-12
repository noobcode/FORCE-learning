function plot_phase_portait_pre_and_post_training(REC, dt, imin, icrit, nt)
    
    pre_time = 1:(imin - round(2000/dt));
    post_time = icrit:(icrit+round(2000/dt));
    
    %% pre training
    for j = 1:1:5
        figure;
        plot(REC(pre_time,j), REC(pre_time,j+5), 'LineWidth', 0.25)
        xlabel('v')
        ylabel('u') 
        title(strcat("Pre Learning - Phase Portrait ", num2str(j)))
    end
    
    %% post training
    for j = 1:1:5
        figure;
        plot(REC(post_time,j), REC(post_time,j+5), 'LineWidth', 0.25) 
        xlabel('v')
        ylabel('u') 
        title(strcat("Post Learning - Phase Portrait ", num2str(j)))
    end
    
    %% all in one figure
    %figure;
    %for i=1:1:5
    %    subplot(5,2,i)
    %    plot(REC(pre_time,j), REC(pre_time,i+5), 'LineWidth', 0.25)
    %    title(strcat("Pre Learning - Phase Portrait ", num2str(j)))
    %end
    
    %for i=1:1:5
        
    %end
    
    
end