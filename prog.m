s = load('er_plots.mat'); er_plots = s.er_plots;
for i = 1:length(er_plots)
    h1 = subplot(1,2,1);
    etalon = load('detector_filtered_field.mat');
    etalon = etalon.detector_filtered_field;
    y2 = er_plots(i).CF; plot(h1,y2,er_plots(i).detect(y2), 'xr');  hold on;
    z2 = er_plots(i).CF_et; plot(h1,z2,etalon(z2), 'xg'); hold on;
    plot(h1,er_plots(i).detect); title(num2str(er_plots(i).number)); hold on;
    plot(h1,etalon); hold off;
    h2 = subplot(1,2,2);
    ERR = er_plots(i).error;
    str = {['Error: ', num2str(ERR)]; ['Ythresh: ', num2str(er_plots(i).k_amp)]; ...
        ['Fluc ne_0: ', num2str(er_plots(i).parameters.p_FlucNe0)]; ...
        ['Fluc I_{cen}: ', num2str(er_plots(i).parameters.p_FlucIcen)]; ...
        ['Fluc J_{cen}: ', num2str(er_plots(i).parameters.p_FlucJcen)]};
    h2 = text(0.5,0.5,str);  
    set(h2, 'Color', 'black', 'HorizontalAlignment', 'center');
    set(gca,'Visible', 'off');
    pause;
    set(h2, 'Visible', 'off');
end;
