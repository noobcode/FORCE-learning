function plot_membrane_pre_and_post_learning(REC, dt, nt, vpeak, vreset, imin, icrit)
    % "normalize" membrane potential,
    % add 'j' to the membrane potential of neuron j so that the curves do not
    % overlap
    pre_time = 1:imin;
    post_time = icrit:nt;
    %% membrane potential
    figure;
    for j = 1:1:5
        plot(pre_time*dt/1000, REC(pre_time,j)/(vpeak-vreset)+1.5*j), hold on 
    end
    %xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    xlabel('Time (s)')
    ylabel('v(t)') 
    set(gca,'YTick',[]);
    title('Pre-Learning')
    
    % post learning
    figure;
    for j = 1:1:5
        plot(post_time*dt/1000, REC(post_time,j)/(vpeak-vreset)+2*j), hold on 
    end
    %xlim([icrit,nt]*dt/1000) % last 5 seconds
    xlabel('Time (s)')
    ylabel('v(t)')
    set(gca,'YTick',[]);
    title('Post Learning')
    
    %% try subplot - membrane
    figure;
    for j = 1:1:5
        subplot(5,1,j)
        plot(pre_time*dt/1000, REC(pre_time,j))
        if j < 5
            set(gca,'XTick',[]);
        end
        set(gca,'YTick',[]);
        ylabel(num2str(j))
    end
    xlabel('Time (s)')
    title('Membrane Potential - Pre Learning')
    %xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    
    figure;
    for j = 1:1:5
        subplot(5,1,j)
        plot(post_time*dt/1000, REC(post_time,j))
        if j < 5
            set(gca,'XTick',[]);
        end
        set(gca,'YTick',[]);
        ylabel(num2str(j))
    end
    xlabel('Time (s)')
    title('Membrane Potential - Post Learning')
    %xlim([icrit,nt]*dt/1000) % from 0 to 5 seconds
    
    %% adaption variable
    
    figure;
    for j = 6:1:10
        plot(pre_time*dt/1000, REC(pre_time,j)), hold on 
    end
    %xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    xlabel('Time (s)')
    ylabel('u(t)')
    title('Pre-Learning')
    
    % post learning
    figure;
    for j = 6:1:10
        plot(post_time*dt/1000, REC(post_time,j)), hold on 
    end
    %xlim([icrit,nt]*dt/1000) % last 5 seconds
    xlabel('Time (s)')
    ylabel('u(t)')
    title('Post Learning')
    
    %% try subplot - adaption
    % pre training
    figure;
    for j = 6:1:10
        subplot(5,1,j-5)
        plot(pre_time*dt/1000, REC(pre_time,j))
        if j < 10
            set(gca,'XTick',[]);
        end
        set(gca,'YTick',[]);
        ylabel(num2str(j-5))
    end
    xlabel('Time (s)')
    title('Adaption Variable - Pre Learning')
    %xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    
    % post training
    figure;
    for j = 6:1:10
        subplot(5,1,j-5)
        plot(post_time*dt/1000, REC(post_time,j))
        if j < 10
            set(gca,'XTick',[]);
        end
        set(gca,'YTick',[]);
        ylabel(num2str(j-5))
    end
    xlabel('Time (s)')
    title('Adaption Variable - Post Learning')
    %xlim([icrit,nt]*dt/1000) 
end