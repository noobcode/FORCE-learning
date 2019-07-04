function orderspikes(N, dt, tspike, startpoint, endpoint)

    %startpoint = 2.418*1000;  % in ms
    %endpoint = 2.4325*1000;
    
    tspike2 = tspike( find( (tspike(:,2) > startpoint) & (tspike(:,2) < endpoint) ),:);
    %tspike2 = tspike( find( (tspike(:,2) < endpoint) ),:);
    
    [x,its] = sort(tspike2(:,2),1);
    clear x;
    itsnum = unique(tspike2(its,1),'stable');
    X=[1:N]';
    Y=find(~ismember(X,itsnum));
    Z=[itsnum;X(Y)];
    M = containers.Map(Z,X);
    
    plotM = [];
    for i=1:length(tspike)
        plotM = [plotM; M(tspike(i,1))];
    end
 
    maxt = max(tspike(:,2)/1000); 
    
    figure;
    plot(tspike(:,2)/1000, plotM,'.');
    xlim([0,maxt]);
end
