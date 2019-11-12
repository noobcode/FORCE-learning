function plot_eigenvalues_pre_and_post_learning(OMEGA, E, BPhi)
    figure;
    Z = eig(OMEGA + E*BPhi'); % eigenvalues after learning 
    Z2 = eig(OMEGA); % eigenvalues before learning 

    plot(Z2,'r.'), hold on 
    plot(Z,'k.') 
    legend('Pre-Learning','Post-Learning')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title('Eigenvalues')
end