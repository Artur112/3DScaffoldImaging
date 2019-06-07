function [th,alive_bar] = new_alive_th(th,data,thicknesses,XLim,t, alive_bar)
    th = round(th,1);
    for week = 1:length(data)
       numalive_lim(week) = sum(data{week}{1}.WHratios > th);
    end
    YLim2 = max(numalive_lim)*1.1;
    for week = 1:length(data)
        alive = [];
        for m = 1:length(data{week})
            alive(m) = sum(data{week}{m}.WHratios > th);
        end  
        ff2 = subplot(3,3,3+week);
        ff2.Position = ff2.Position + [0, 0.02,0,0];
        distance = (1:length(alive)).*str2double(thicknesses{week});
        bar(distance,alive);
        alive_bar(week) = gca;
        if(week == 2)
            title('Number of Alive Cells');
        end
        if(week == 1)
         ylabel('Number of alive cells');
        end
        ylim([0 YLim2]);
        xlim([0 XLim]);
    end
    t.String = {'Width to Height Ratio Threshold:',num2str(th)};
end