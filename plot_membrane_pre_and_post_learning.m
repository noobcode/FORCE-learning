function plot_membrane_pre_and_post_learning(REC, dt, nt, vpeak, vreset, imin, icrit)
    % "normalize" membrane potential,
    % add 'j' to the membrane potential of neuron j so that the curves do not
    % overlap
    figure;
    for j = 1:1:5
        plot((1:1:nt)*dt/1000, REC(1:1:nt,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    xlabel('Time (s)')
    ylabel('Neuron Index') 
    title('Pre-Learning')
    
    % post learning
    figure;
    for j = 1:1:5
        plot((1:1:nt)*dt/1000, REC(1:1:nt,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([icrit,nt]*dt/1000) % last 5 seconds
    xlabel('Time (s)')
    ylabel('Neuron Index') 
    title('Post Learning')
end