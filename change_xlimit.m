function x_limit = change_xlimit(x_limit,cell_bar,alive_bar,scaff_bar,numweeks,thicknesses,ax2,t3,scaffperc)
    x_limit = floor(x_limit);
    % Change x limits of all the plots
    for week = 1:numweeks
        set(cell_bar(week),'XLim',[0 x_limit]);
        set(alive_bar(week),'XLim',[0 x_limit]);
        set(scaff_bar(week),'XLim',[0 x_limit]);
        set(ax2(week),'XLim',[0, round(x_limit/str2double(thicknesses{week}))]...
            ,'XTick',0:round(x_limit/str2double(thicknesses{week})/5):round(x_limit/str2double(thicknesses{week})));
    end
    % Make the ylimits the same for all the scaffold plots
    for week = 1:numweeks
        lim = round(x_limit/str2double(thicknesses{week}));
        if(lim < length(scaffperc{week}))
            scaff_yLimits(week) = max(scaffperc{week}(1:round(x_limit/str2double(thicknesses{week}))));  %# Get the range of the y axis
        else
            scaff_yLimits(week) = max(scaffperc{week});
        end
    end
    new_scaff_lim = max(scaff_yLimits)*1.1;
    for week = 1:numweeks
        set(scaff_bar(week),'YLim',[0 new_scaff_lim]);
    end
    t3.String = num2str(x_limit) + " \bf\mu\itm";
end