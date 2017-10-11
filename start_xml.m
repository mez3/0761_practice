start = 1;
stop = 99;
step = 1;
xml = zeros(1,stop-start);
xml_min = 1000000;
%j = 1;
%k = 0;
%xml_min(1) = 1000000;
%xml_min(2) = 1000000;
% while k < 5
%     j = j + 1;
%     for i = start : step : stop
%         xml(i) = tof_mat(i,  2 * j + 1);
%         %pr = i * 100 / (stop - start);
%         text = [num2str(i), '% : ', num2str(xml(i))];
%         disp(text);
% %         if xml(i) < xml_min
% %             xml_min(j) = xml(i);
% %             i_min(j) = i;
% %         end;
%     end;
%     [xml_min(j), i_min(j)] = min(xml);
%     if xml_min(j) > xml_min(j-1)
%         k = k + 1;
%     end;
% end;
for i = start : step : stop
        xml(i) = tof_mat(i);
        %pr = i * 100 / (stop - start);
         text = [num2str(i), ' : ', num2str(xml(i))];
         disp(text);
          if xml(i) < xml_min
            xml_min = xml(i);
            i_min = i;
         end;
end;

%text1 = [num2str(2 * i_absmin + 1), ' : ', num2str(absmin), ' : ', num2str(i_min(i_absmin))];
text = [num2str(i_min), ' ', num2str(xml_min)];
disp(text);
%save('xml.xls', 'xml_min');