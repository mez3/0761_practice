function [k_errors] = tof_mat(k_amp)
    MED=50;% tof_secder median filter windowSz FiltWindSz.MED
    TOFCF_ID_IDX=1;
    visopt = 0;
    delay_cf=1000;
    delay_sd=round(0.5*MED);
    dt=1e-3;
    t_start=0;
    t_end=15;
    shft = 3000;
    etalon = load('detector_filtered_field.mat');
    d=load('det_struct_2.mat');det_struct=d.det_struct;
    n = size(det_struct,2);
    t = t_start:dt:t_end - dt;
    data_et(:,2) = etalon.detector_filtered_field;
    data_et(:,1)=t*1e-9;
    
    
    for j = 2:n
        s = det_struct(j).detect;
        t = t_start:dt:t_end - dt;
        data(:,2)=s; 
        data(:,1)=t*1e-9;
        data2(:,2) = signal_cut(data(:,2));
        data2(:,1) = signal_cut(data(:,1));
       % medfilt1(data2, rang_filt);
        %[id_sd,der2x_f] = tof_secder(data2, shft, MED, delay_sd );
       [id_cf,f_res,y_thresh] = tof_cf(data2,delay_cf,TOFCF_ID_IDX);
       det_struct(j).f_res = f_res;
       det_struct(j).y_thresh = y_thresh;
         data3(:,2) = signal_cut(data_et(:,2));
         data3(:,1) = signal_cut(data_et(:,1));
%         %[id_et_sd,der2x_f_et] = tof_secder(data3, shft, MED, delay_sd );
         [id_et_cf,f_res_et] = tof_cf(data3,delay_cf,TOFCF_ID_IDX);
        %det_struct(j).det_SD_m = id_sd;
        det_struct(j).det_CF_m = id_cf;
        %det_struct(j).det_SD_et = id_et_sd;
      det_struct(j).det_CF_et = id_et_cf;
    end
   % for_approx = det_struct(896).detect;
    save('det_struct_1_m.mat', 'det_struct');
   % save('det_for_approx.mat', 'for_approx');
    index = 100;
    k_errors = visualize();
    
    function [y] = signal_cut(s)
        shift = 1000;
        s_shift = s(shift:end);  
        s_fill = s(shift) * ones(shift-1,1);
        y = [s_fill;s_shift];
    end
    
    function [f_res,f1,f2]=constant_fraction(y,delay)
            f = 0.5;
            EndVal=max(y);
            offset = y(delay)*ones(delay, 1);
            offset2= y(end)*ones(delay,1);
            f1 = vertcat(offset, y, EndVal);
            f2 = f * vertcat(y, offset2,EndVal);
            f_res = f1 - f2;

    end

    function [id_cf,f_res,y] = tof_cf(s,delay,ID_IDX)
        
        Ythreshold=k_amp*0.01*max(s(:, 2)); % 5e-4 % 0.35, 0.3,0.2,0.1,0.08,0.0625,0.05,0.025,0.012

        y = s(:, 2);

        % constant fraction %
        THRESH_VARIANT = 1;
        switch THRESH_VARIANT
            case 1
                y=y-Ythreshold;
                f_res=constant_fraction(y,delay);
            case 2
                y=y-Ythreshold; idx=find(y<0); y(idx)=0;
                f_res=constant_fraction(y,delay);
            case 3
                f_res=constant_fraction(y,delay);
                f_res=f_res-Ythreshold; 
            case 4 % F_RES_AUTO_CORRECTION  
                f_res=constant_fraction(y,delay);
                f_res =f_res - f_res(1)-max(f_res/250); % 5e-4;  
            case 5 % bad method
                [f_res,f1,f2]=constant_fraction(y,delay);
                f2=f2-Ythreshold; idx=find(f2<0); f2(idx)=0;
                f_res = f1 - f2;
        end
        
        % Zerro-crossing %
        t1 = size(f_res, 1);
        k=1;
        for i = 2 + delay : t1
            if f_res(i) * f_res(i - 1) < 0
                id(k) = i - delay;  k = k + 1;
            end
        end

        %  Making id_cf 
            %id_cf=id(k-1)-1;
        k=k-1;
        if ~isempty(id) 
            if ID_IDX > k
                if id(k) > 1
                   id_cf=id(k);
                else
                    id_cf=delay;
                end
            else
                id_cf=id(ID_IDX);
            end
        else
            id_cf=length(y);
        end

        if id_cf >= size(f_res, 1)
           id_cf=length(y); 
        else 
            id_cf=id_cf-1;
        end


    end

    function [id,der2s_f] = tof_secder(s,shft,windowSize,delay)
        SECDER_ID_IDX=1;
        release=1;
        switch release
            case 1
            s=s(shft:end,1:2);
            t=s(:,1);
            y=s(:,2);

            y=[diff(y,2)]; % second derivative
            y= filtfilt(ones(1, windowSize) / windowSize, 1, y); % smoothing

            if 0 
                y=[diff(y,1)]; % Zero-crossing
                y= filtfilt(ones(1, windowSize) / windowSize, 1, y); % smoothing
            end

            idx=find(y>0); y(idx)=0; y=abs(y); % select first pick of choosen sign

            y=y./max(y); % normalization
            der2s_f = [t(1:length(y),1), y];
            id=delay+shft+tof_cf(der2s_f,delay,SECDER_ID_IDX); % bounding to pulse

            if 0    % another not so good id-finder 
                ind=(find(y>0.1*max(y))); 
                id=round((windowSize-1)/2)+shft+ind(1);
            end

            case 2
            s=s(shft:end,1:2);
            t=s(:,1);
            y=s(:,2);

            y=[diff(y,2)]; % second derivative
            y= filtfilt(ones(1, windowSize) / windowSize, 1, y);% smoothing
            % Zero-crossing
                t1 = size(y, 1);
                k=1; idd=[];
                for i = 2+delay:t1
                    if y(i) * y(i - 1) < 0
                        idd(k) = i  - delay; k=k+1;
                    else 
                        %continue
                    end
                end
            % Zero-crossing id-postprocessing    
                % id_cf=id(k-1)-1;
                k=k-1; ID_IDX=1;
                if ~isempty(idd) 
                    if ID_IDX > k
                        if idd(k) > 1
                           id_cs=idd(k);
                        else
                            id_cs=delay;
                        end
                    else
                        id_cs=idd(ID_IDX);
                    end
                else
                    id_cs=length(y);
                end

                if id_cs >= size(y, 1)
                   id_cs=length(y); 
                else 
                    id_cs=id_cs-1;
                end
            %

            if 1
                der2s_f = [t(1:length(y),1), y];
                id=delay+shft+id_cs;    
            else
                der2s_f = [t(1:length(y),1), y];
                id=delay+shft+tof_cf(der2s_f,delay,SECDER_ID_IDX);    
            end

        end

    end
    function [k_er] = visualize()
        %y1 = zeros(size(det_struct, 2) - 1,1);
        y2 = zeros(size(det_struct, 2) - 1,1);
        for i=2:size(det_struct,2)
            %y1(i) = det_struct(i).det_SD_m;
            y2(i) = det_struct(i).det_CF_m;
            %z1(i) = det_struct(i).det_SD_et;
            %z2(i) = det_struct(i).det_CF_et;
        end
        %m1 = mean(y1)-2000;
        %for_gr = find(y1 > m1);
       % m_y1 = mean(y1(1:50));
        m_y2 = mean(y2(2:5));
        er_y2_1 = find(y2 < (m_y2-500));
        er_y2_2 = find(y2 > (m_y2+500));
      
   
        errors2 = vertcat(er_y2_1, er_y2_2);
        errors2 = errors2';
        n_er = length(errors2);
        er_plots(n_er,1) = struct('detect','','number','','CF','','error','','k_amp','');
        k = 1;
        for i = errors2
            if i > 1
                er_plots(k).detect = det_struct(i).detect;
                er_plots(k).CF = det_struct(i).det_CF_m;
                er_plots(k).number = i;
                er_plots(k).error = abs(y2(i)-m_y2);
                er_plots(k).k_amp = k_amp;
                er_plots(k).parameters = det_struct(i).parameters;
                er_plots(k).CF_et = det_struct(i).det_CF_et;
                er_plots(k).f_res = det_struct(i).f_res;
                er_plots(k).y_thresh = det_struct(i).y_thresh;
               
                k = k + 1;
            end;
        end;
        
        save('er_plots.mat', 'er_plots');
        k_er = n_er*100/length(y2);
%         for i = idx
%        
%             h1=subplot(2,2,1);
%  
%                 y1(i) = det_struct(i).det_SD_m; plot(h1,y1(i),det_struct(i).detect(y1(i)),'or'); hold on;
%                 y2(i) = det_struct(i).det_CF_m; plot(h1,y2(i),det_struct(i).detect(y2(i)), 'xg'); hold on;
%                 z1(i) = det_struct(i).det_SD_et; plot(h1,z1(i),et_signal(z1(i)),'or'); hold on;
%                 z2(i) = det_struct(i).det_CF_et; plot(h1,z2(i),et_signal(z2(i)), 'xg'); hold on;
%                 plot(h1,et_signal, 'b'); hold on;
%                 plot(h1,det_struct(i).detect);title([num2str(i),' : ',num2str(abs(m_y1-y1(i))), ' : ', num2str(abs(m_y2-y2(i)))]); hold off; 
%             h2=subplot(2,2,2);
%                p_plotFlucDropouts(h2,i);
%             pause();
%         end
if(visopt)        
        figure;
        
        % plot(mt, m_y1, '-'); hold on;
         plot( m_y2 * ones(1, length(y2)), '-'); hold on;
        %plot(y1, '.r'); hold on;
        plot(y2, '.g'); hold on;
        %plot(z1, '.b'); hold on;
       % plot(z2, '.c'); hold off;
end;
        
    end

    function [xPolygon,yPolygon]=p_calcPolygonBound(i)
        p_Imax=det_struct(i).parameters.p_Imax;
        p_Jmax=det_struct(i).parameters.p_Jmax;
        yPolygon=[0,p_Imax,p_Imax,0,0]; xPolygon=[0,0,p_Jmax,p_Jmax,0];
    end
            
    function [xCircle1,yCircle1]=p_calcFlucBound(i)
        p_FlucIcen=det_struct(i).parameters.p_FlucIcen;
        p_FlucJcen=det_struct(i).parameters.p_FlucJcen;
        p_FlucRadius=det_struct(i).parameters.p_FlucRadius;

        ii=0;
        for phi=0:10:360
            ii=ii+1;
            xCircle1(ii)=p_FlucJcen+p_FlucRadius*cos(phi/180*pi);
            yCircle1(ii)=p_FlucIcen+p_FlucRadius*sin(phi/180*pi);
        end
    end
    
    function p_plotFlucDropouts(h,i)
       [xPolygon,yPolygon]=p_calcPolygonBound(2);
        plot(h,xPolygon,yPolygon,'b'); hold on
       [xCircle1,yCircle1]=p_calcFlucBound(i);
        plot(h,xCircle1,yCircle1,'-r');
        daspect([1 1 1]);hold off
    end

    

end
