s = load('er_plots.mat'); er_plots = s.er_plots; 
k = 1;
for i = 1:length(er_plots)
    if i > 1
        if (er_plots(i).number - er_plots(i-1).number) == 1
            count_group = count_group + 1;
        else 
            len_group = count_group;
            count_group = 1;
            if len_group > max_group
                max_group = len_group;
                number_mg = er_plots(i-1).number - len_group;
            end;
        end;
    else
        count_group = 1;
        len_group = 1;
        max_group = 1;
        number_mg = 1;
    end;
    
    h1 = subplot(1,2,1);
    etalon = load('detector_filtered_field.mat');
    etalon = etalon.detector_filtered_field;
    y2 = er_plots(i).CF; plot(h1,y2,er_plots(i).detect(y2), 'xr');  hold on;
    z2 = er_plots(i).CF_et; plot(h1,z2,etalon(z2), 'xg'); hold on;
   % yv = er_plots(i).CF_et * ones(1,5000); plot(h1, yv,'g'); hold on;
    plot(h1,er_plots(i).detect); title(num2str(er_plots(i).number)); hold on;
    plot(h1,etalon); hold on;
    plot(h1, er_plots(i).y_thresh); hold off;
    h2 = subplot(1,2,2);
    ERR = er_plots(i).error;
    str = {[num2str(i),' / ', num2str(length(er_plots)) ];['Error: ', num2str(ERR)]; ...
        ['Ythresh: ', num2str(er_plots(i).k_amp)]; ...
        ['Fluc ne_0: ', num2str(er_plots(i).parameters.p_FlucNe0)]; ...
        ['Fluc I_{cen}: ', num2str(er_plots(i).parameters.p_FlucIcen)]; ...
        ['Fluc J_{cen}: ', num2str(er_plots(i).parameters.p_FlucJcen)]; ...
        ['Number in current group: ', num2str(count_group)]; ['Length of previous group: ', num2str(len_group)]; ...
        ['Max length of group: ', num2str(max_group)]; ['Place of max group: ', num2str(number_mg)]}; 
    h2 = text(0.5,0.5,str);  
    set(h2, 'Color', 'black', 'HorizontalAlignment', 'center', 'FontSize', 18);
    set(gca,'Visible', 'off');
    pause;
    set(h2, 'Visible', 'off');
end;
 