function call_orderspikes(tspike, f, t_start, t_end)
    period_length = 1000/f; % period length, ms
    
    start_points = t_start:period_length:t_end;
    end_points = start_points + period_length;
    
    for i=1:length(start_points)
        orderspikes(2000, 0.04, tspike, start_points(i), end_points(i))
    end

end