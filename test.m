%%-------------------------------------------------------------------------
% Created by    : Suba Siva Chandran
% Last edited on: 19/02/2021
% Description   : To find the degree of cure from the experimental curves,
%                 this script has been created where it can find the initial
%                 and final peak temperatures. With those temperatures the
%                 total area can be calculated. Also the user can enter any
%                 temperature so as to find the degree of cure at that point.
%                 The calculations has been made according to the standard DIN
%                 51007:2019-04 (General Principles of Differential thermal 
%                 analysis (DTA) and differential scanning calorimetry
%                 (DSC)).
%%-------------------------------------------------------------------------
clear all;
clc;
Location = uigetdir;
% Identify where to search for files
% Store the name of all .xls files as a vector D
% D = dir([Location, '\*.xlsx']);
source_files = dir(fullfile(Location, '*.xlsx'));
dstfile = 'results.xlsx';
len= length(source_files);
looptables = cell(len, 1);
for i=1:len
    clear final_peak
    data= xlsread(source_files(i).name);
% Extract the file names
% filenames = {D(:).name}.';
% data = cell(length(D),1);
% for ii = length(D):-1:1 
%       % Create the full file name and partial filename
%       fullname = [Location filesep D(ii).name];
%       % Read in the data
%       data{ii} = xlsread(fullname);         
% end


%Loading the data from spreadsheet
% data = readtable('testdata1.xlsx');
% x = data.Temp;
% y = data.DSC; 
x= data(:,1);
y= data(:,2);
[~,max_idx] = max(y);
peak_temp = data(max_idx,:);
Max_temperature = peak_temp(:,1);
Max_DSC = peak_temp(:,2);
% [~,max_idx] = max(data.DSC);
% peak_temp = data(max_idx,:);
% Max_temperature = peak_temp.Temp;
% Max_DSC = peak_temp.DSC;
m = numel(x);

% Find slope to calculate initial and final temperature values
slope_new = [];
[gy] = gradient(y);
[gx] = gradient(x);
for k=1:m
   if x(k) >= 80 & x(k) <= 200 %to calculate only between temperatures 80?C and 200?C (user can vary the range to calculate the slope)
       slope_temp = [(gy(k) ./ gx(k)),x(k),y(k)];
       slope_new = [slope_new;slope_temp];
   end
end

%Finding the Initial peak temperature value using the slope
n = numel(slope_new(:,1));
for l =1:n
if slope_new(l) >= 0.01
    initial_peak = [slope_new(l,2),slope_new(l,3)];
    break
end
end

%Finding the Final peak temperature value using the slope
p = numel(slope_new(:,2));
for o =1:p
               if  slope_new(o,2) > Max_temperature
                    while slope_new(o,1) >= 0
                    final_peak = [slope_new(o,2),slope_new(o,3)];
                    break
                    end
               end
if exist ('final_peak') == 1
break
end
end

%Calculated Initial and Final peak temperature and DSC values
initial_peaktemp = initial_peak(:,1);
initial_peakDSC = initial_peak(:,2);
final_peaktemp = final_peak(:,1);
final_peakDSC = final_peak(:,2);

%%%To generate a linear line between initial and final peak values
x_index_points = x>= initial_peaktemp & x<=final_peaktemp;
x_index_linear = x(x_index_points);
x_index_linear_count = numel(x_index_linear);
y_index_linear = linspace(initial_peakDSC,final_peakDSC,x_index_linear_count);
% plot(x,y,x_index_linear,y_index_linear,'--');

%%To find x where y = 0 to calculate the area
y_0 = y*0; %To plot the horizontal line along 0

%%%%Find x-intercept at y=0
xc = find((circshift(slope_new(:,3), [0 1])) > 0);
length_xc = length(xc);
xc_y_begin = slope_new(xc(1)-1,3);
xc_y_end = slope_new(xc(length(xc))+1,3);
xc_y = [xc_y_begin;slope_new(xc,3);xc_y_end];
length_xc_y = length(xc_y);
xc_x_begin = slope_new(xc(1)-1,2);
xc_x_end = slope_new(xc(length(xc))+1,2);
xc_x = [xc_x_begin;slope_new(xc,2);xc_x_end];
length_xc_x = length(xc_x);
% xc = find((circshift(y, [0 1])) > 0); % Define ?y? Zero-Crossings
idx_begin = xc(1)-1:xc(1);% Index Range Including Zero-Crossing
idx_end = xc(length_xc):xc(length_xc)+1;
x_begin = interp1(xc_y(1:2), xc_x(1:2), 0);        % Find ?x? At Zero Crossing
x_end = interp1(xc_y(length_xc_y-1:length_xc_y), xc_x(length_xc_x-1:length_xc_x), 0);
% x_begin = interp1(y(xc(1)-1:xc(1)), x(xc(1)-1:xc(1)), 0);        % Find ?x? At Zero Crossing
% x_end = interp1(y(xc(length_xc):xc(length_xc)+1), x(xc(length_xc):xc(length_xc)+1), 0);

index = find(slope_new(:,3)>=0);
x_point_index = slope_new(index,:);
x_point = x_point_index(1,2);
x_point_DSC = x_point_index(1,3);
index_2 = find( slope_new(:,2) > Max_temperature & slope_new(:,3)<= 0);
x_point2_index = slope_new(index_2,:);
x_point2 = x_point2_index(1,2);
x_point2_DSC = x_point2_index(1,3);

%Calculating the area of the upper and lower part
gt0_x = x>= x_begin & x<= x_end & y>=0;
xu = x(gt0_x);
xu = [x_begin;xu;x_end];
yu = y(gt0_x);
yu = [0;yu;0];
lt0_x = x>= x_end &  x <= final_peaktemp & y<=0 ;
lt0_x_2 = x>= initial_peaktemp & x <= x_begin & y<=0;
xl1 = x(lt0_x);
xl1 = [x_end;xl1];
yl1 = y(lt0_x);
yl1 = [0;yl1];
xl2 = x(lt0_x_2);
xl2 = [xl2;x_begin];
yl2 = y(lt0_x_2);
yl2 = [yl2;0];
posArea = trapz(xu,yu);
% posArea1 = sum(cumtrapz(xu,yu));
x_coordinates = [initial_peaktemp initial_peaktemp final_peaktemp final_peaktemp initial_peaktemp];
y_coordinates = [initial_peakDSC 0 0 final_peakDSC initial_peakDSC];
negArea1 = trapz(xl1,yl1);
negArea2 = trapz(xl2,yl2);
negArea3 = trapz(x_coordinates,y_coordinates);
negArea_final = abs(negArea1 + negArea2 + negArea3);
Total_area = posArea + negArea_final
% Max_temperature = peak_temp.Temp
Max_temperature = peak_temp(:,1)
initial_peaktemp = initial_peak(:,1)
final_peaktemp = final_peak(:,1)

%%%Plot to find the total area
figure
plot(x,y,Max_temperature,Max_DSC,'ro',initial_peaktemp,initial_peakDSC,'rx',final_peaktemp,final_peakDSC,'rx',x,y_0,'-',x_begin,0,'bx',x_end,0,'bx',x_index_linear,y_index_linear,'-');
hold on
text(Max_temperature,Max_DSC,sprintf('Max.temp=%f',Max_temperature));
text(initial_peaktemp,initial_peakDSC,sprintf('Ini.peak=%f',initial_peaktemp));
text(final_peaktemp,final_peakDSC,sprintf('Fin.peak=%f',final_peaktemp));
title('Area under the curve calculation')
xlabel('Temperature (\circC)')
ylabel('DSC (mW/mg)')
a1 = fill(xu,yu,'r');
a2 = area(xl1,yl1);
a3 = area(xl2,yl2);
a3.FaceColor = [0 0.5 1];
% a4 = fill(x_coordinates,y_coordinates,'y');
A1str = sprintf('Area 1 = %6.3f', posArea);
A2str = sprintf('Area 2 = %6.3f', negArea1);
A3str = sprintf('Area 3 = %6.3f', negArea2);
Ini.peak = sprintf('Ini.peak.temp = %6.3f',initial_peaktemp);
% A4str = sprintf('Area 4 = %6.3f', negArea3);
legend([a1 a2 a3],A1str,A2str,A3str)
hold off

%%Dialog box for user to enter any temperature to find the degree of cure
prompt = {'Enter the temperature value to find degree of cure (in celsius)'};
dlgtitle = 'Temperature Value';
definput = {' '};
dims = [1 70];
% opts.Interpreter = 'tex';
anytemperature = inputdlg(prompt,dlgtitle,dims,definput);
Beliebige_temp = str2num(anytemperature{1});

%%%%Find y-intercept at x= Beliebige_temp
%%%Upper y-intercept
yc = find((circshift(x, [0 1])) > Beliebige_temp);
length_yc = length(yc);
yc_x_begin = x(yc(1)-1,1);
yc_x_end = x(yc(1));
yc_y_begin = y(yc(1)-1,1);
yc_y_end = y(yc(1));
yc_x = [yc_x_begin;yc_x_end];
yc_y = [yc_y_begin;yc_y_end];
y_intercept_bel_temp_up = interp1(yc_x(1:2),yc_y(1:2), Beliebige_temp);

%%%Lower y-intercept
xy_bel_temp_points = [x_index_linear,y_index_linear'];
loweryc = find((circshift(xy_bel_temp_points(:,1), [0 1])) > Beliebige_temp);
yc_x_begin_low = xy_bel_temp_points(loweryc(1)-1,1);
yc_x_end_low = xy_bel_temp_points(loweryc(1),1);
yc_y_begin_low = xy_bel_temp_points(loweryc(1)-1,2);
yc_y_end_low = xy_bel_temp_points(loweryc(1),2);
yc_x_low = [yc_x_begin_low;yc_x_end_low];
yc_y_low = [yc_y_begin_low;yc_y_end_low];
y_intercept_bel_temp_low = interp1(yc_x_low(1:2),yc_y_low(1:2), Beliebige_temp);

%Finding area for beliebige temperature for y>0 between x_begin and x_end
if Beliebige_temp >= x_begin & Beliebige_temp <= x_end
gt0_x_1 = x>= x_begin & x<= Beliebige_temp &  y>=0 ;
xu1 = x(gt0_x_1);
xu1 = [x_begin;xu1;Beliebige_temp;Beliebige_temp];
yu1 = y(gt0_x_1);
yu1 = [0;yu1;y_intercept_bel_temp_up;0];
Bel_temp_pos_area = trapz(xu1,yu1);

%Finding area for beliebige temperature for y<0 between x_begin and x_end
x_bel_coord_low = [initial_peaktemp initial_peaktemp Beliebige_temp Beliebige_temp initial_peaktemp];
y_bel_coord_low = [initial_peakDSC 0 0 y_intercept_bel_temp_low initial_peakDSC];
xy_bel_low_area = trapz(x_bel_coord_low,y_bel_coord_low);
Bel_temp_neg_area = negArea2 + xy_bel_low_area;
Bel_temp_area = Bel_temp_pos_area + Bel_temp_neg_area;

%%%Degree of cure
Deg_of_cure = Bel_temp_area / Total_area *100

%%%Plot for beliebige temperature
figure
plot(x,y,Max_temperature,Max_DSC,'ro',initial_peaktemp,initial_peakDSC,'rx',final_peaktemp,final_peakDSC,'rx',x,y_0,'-',x_begin,0,'bx',x_end,0,'bx',x_index_linear,y_index_linear,'-',Beliebige_temp,0,'gx',Beliebige_temp,y_intercept_bel_temp_up,'gx',Beliebige_temp,y_intercept_bel_temp_low,'gx');
y_limit=get(gca,'ylim');
hold on
text(Max_temperature,Max_DSC,sprintf('Max.temp=%f',Max_temperature));
text(initial_peaktemp,initial_peakDSC,sprintf('Ini.peak=%f',initial_peaktemp));
text(final_peaktemp,final_peakDSC,sprintf('Fin.peak=%f',final_peaktemp));
title('Estimation curve for finding the curing rate')
xlabel('Temperature (\circC)')
ylabel('DSC (mW/mg)')
plot([Beliebige_temp Beliebige_temp], y_limit);
a1_Bel_area = fill(xu1,yu1,'r');
A4str = sprintf('Degree of cure = %6.3f', Deg_of_cure);
legend([a1_Bel_area],A4str)
hold off

%Finding area for beliebige temperature for x < x_begin(x-intercept at y = 0)
elseif Beliebige_temp <= x_begin
xy_bel_temp_points = [x_index_linear,y_index_linear'];
loweryc = find((circshift(xy_bel_temp_points(:,1), [0 1])) > Beliebige_temp);
yc_x_begin_low = xy_bel_temp_points(loweryc(1)-1,1);
yc_x_end_low = xy_bel_temp_points(loweryc(1),1);
yc_y_begin_low = xy_bel_temp_points(loweryc(1)-1,2);
yc_y_end_low = xy_bel_temp_points(loweryc(1),2);
yc_x_low = [yc_x_begin_low;yc_x_end_low];
yc_y_low = [yc_y_begin_low;yc_y_end_low];
y_intercept_bel_temp_low = interp1(yc_x_low(1:2),yc_y_low(1:2), Beliebige_temp);

x_bel_coord_low = [initial_peaktemp initial_peaktemp Beliebige_temp Beliebige_temp initial_peaktemp];
y_bel_coord_low = [initial_peakDSC 0 0 y_intercept_bel_temp_low initial_peakDSC];
xy_bel_low_area = trapz(x_bel_coord_low,y_bel_coord_low);
xbegin_lt0_x = x>= initial_peaktemp & x <= Beliebige_temp;
negArea_bel_2_xcoord = x(xbegin_lt0_x);
% negArea_bel_2_xcoord_plot_linear = find((circshift(xy_bel_temp_points(:,1), [0 1])) < Beliebige_temp);
% negArea_bel_2_xcoord_plot_1 = xy_bel_temp_points(negArea_bel_2_xcoord_plot_linear,1);
% negArea_bel_2_xcoord_plot = [negArea_bel_2_xcoord;Beliebige_temp;Beliebige_temp;negArea_bel_2_xcoord_plot_1];
negArea_bel_2_xcoord = [negArea_bel_2_xcoord;Beliebige_temp;Beliebige_temp];
negArea_bel_2_ycoord = y(xbegin_lt0_x);
% negArea_bel_2_ycoord_plot_1 = xy_bel_temp_points(negArea_bel_2_xcoord_plot_linear,2);
% negArea_bel_2_ycoord_plot = [negArea_bel_2_ycoord;y_intercept_bel_temp_up;y_intercept_bel_temp_low;negArea_bel_2_ycoord_plot_1];
negArea_bel_2_ycoord = [negArea_bel_2_ycoord;y_intercept_bel_temp_up;0];
negArea_bel_2 = trapz(negArea_bel_2_xcoord,negArea_bel_2_ycoord);
Bel_temp_area = negArea_bel_2 + xy_bel_low_area; 
Deg_of_cure = Bel_temp_area / Total_area *100

%%%Plot for beliebige temperature when x < x_begin
figure
plot(x,y,Max_temperature,Max_DSC,'ro',initial_peaktemp,initial_peakDSC,'rx',final_peaktemp,final_peakDSC,'rx',x,y_0,'-',x_begin,0,'bx',x_end,0,'bx',x_index_linear,y_index_linear,'-',Beliebige_temp,0,'gx',Beliebige_temp,y_intercept_bel_temp_up,'gx',Beliebige_temp,y_intercept_bel_temp_low,'gx');
y_limit=get(gca,'ylim');
hold on
text(Max_temperature,Max_DSC,sprintf('Max.temp=%f',Max_temperature));
text(initial_peaktemp,initial_peakDSC,sprintf('Ini.peak=%f',initial_peaktemp));
text(final_peaktemp,final_peakDSC,sprintf('Fin.peak=%f',final_peaktemp));
title('Estimation curve for finding the curing rate')
xlabel('Temperature (\circC)')
ylabel('DSC (mW/mg)')
plot([Beliebige_temp Beliebige_temp], y_limit);
xl3_coordinates = [initial_peaktemp Beliebige_temp Beliebige_temp initial_peaktemp];
yl3_coordinates = [initial_peakDSC y_intercept_bel_temp_up y_intercept_bel_temp_low initial_peakDSC];
% a1_Bel_area = fill(xl3_coordinates,yl3_coordinates,'r');
% fillStart = find(x>=initial_peaktemp,1);
% fillEnd = find(x>=Beliebige_temp,1);
% area(negArea_bel_2_xcoord_plot,negArea_bel_2_ycoord_plot,negArea_bel_2_ycoord_plot_1,'FaceColor', [0.5 0.5 0.5]);
% area(x(fillStart:fillEnd),y(fillStart:fillEnd),initial_peakDSC,'FaceColor', [0.5 0.5 0.5]);
a4 = area(negArea_bel_2_xcoord,negArea_bel_2_ycoord,'FaceColor', [1 0 0]);
% a4 = area(xl3_coordinates,yl3_coordinates);
A5str = sprintf('Degree of cure = %6.3f', Deg_of_cure);
legend([a4],A5str)
hold off

%Finding area for beliebige temperature for x > x_end(x-intercept at y = 0)
elseif Beliebige_temp >= x_end
xy_bel_temp_points = [x_index_linear,y_index_linear'];
loweryc = find((circshift(xy_bel_temp_points(:,1), [0 1])) > Beliebige_temp);
yc_x_begin_low = xy_bel_temp_points(loweryc(1)-1,1);
yc_x_end_low = xy_bel_temp_points(loweryc(1),1);
yc_y_begin_low = xy_bel_temp_points(loweryc(1)-1,2);
yc_y_end_low = xy_bel_temp_points(loweryc(1),2);
yc_x_low = [yc_x_begin_low;yc_x_end_low];
yc_y_low = [yc_y_begin_low;yc_y_end_low];
y_intercept_bel_temp_low = interp1(yc_x_low(1:2),yc_y_low(1:2), Beliebige_temp);

% x_bel_coord_low = [initial_peaktemp x_begin Beliebige_temp Beliebige_temp initial_peaktemp];
% y_bel_coord_low = [initial_peakDSC 0 0 y_intercept_bel_temp_low initial_peakDSC];
x_bel_coord_range = x >= initial_peaktemp & x<= x_begin;
x_bel_coord_low = x(x_bel_coord_range);
x_bel_coord_low = [x_bel_coord_low;x_begin;Beliebige_temp;Beliebige_temp;initial_peaktemp];
y_bel_coord_low = y(x_bel_coord_range);
y_bel_coord_low = [y_bel_coord_low;0;0;y_intercept_bel_temp_low;initial_peakDSC];
xy_bel_low_area = trapz(x_bel_coord_low,y_bel_coord_low);
xbegin_lt0_x = x<= Beliebige_temp & x >= x_end;
negArea_bel_2_xcoord = x(xbegin_lt0_x);
negArea_bel_2_xcoord = [x_end;negArea_bel_2_xcoord;Beliebige_temp;Beliebige_temp];
negArea_bel_2_ycoord = y(xbegin_lt0_x);
negArea_bel_2_ycoord = [0;negArea_bel_2_ycoord;y_intercept_bel_temp_up;0];
negArea_bel_2 = trapz(negArea_bel_2_xcoord,negArea_bel_2_ycoord);
Bel_temp_area = negArea_bel_2 + xy_bel_low_area; 
Final_areaundercurve = posArea + Bel_temp_area;
Deg_of_cure = Final_areaundercurve / Total_area *100

%%%Plot for beliebige temperature when x > x_end
figure
plot(x,y,Max_temperature,Max_DSC,'ro',initial_peaktemp,initial_peakDSC,'rx',final_peaktemp,final_peakDSC,'rx',x,y_0,'-',x_begin,0,'bx',x_end,0,'bx',x_index_linear,y_index_linear,'-',Beliebige_temp,0,'gx',Beliebige_temp,y_intercept_bel_temp_up,'gx',Beliebige_temp,y_intercept_bel_temp_low,'gx');
y_limit=get(gca,'ylim');
hold on
text(Max_temperature,Max_DSC,sprintf('Max.temp=%f',Max_temperature));
text(initial_peaktemp,initial_peakDSC,sprintf('Ini.peak=%f',initial_peaktemp));
text(final_peaktemp,final_peakDSC,sprintf('Fin.peak=%f',final_peaktemp));
title('Estimation curve for finding the curing rate')
xlabel('Temperature (\circC)')
ylabel('DSC (mW/mg)')
plot([Beliebige_temp Beliebige_temp], y_limit);
xl3_coordinates = [initial_peaktemp Beliebige_temp Beliebige_temp initial_peaktemp];
yl3_coordinates = [initial_peakDSC y_intercept_bel_temp_up y_intercept_bel_temp_low initial_peakDSC];
% a1_Bel_area = fill(xl3_coordinates,yl3_coordinates,'r');
% fillStart = find(x>=initial_peaktemp,1);
% fillEnd = find(x>=Beliebige_temp,1);
% area(x(fillStart:fillEnd),y(fillStart:fillEnd),'FaceColor', [0.5 0.5 0.5]);
a1 = fill(xu,yu,'r');
% a4 = area(negArea_bel_2_xcoord,negArea_bel_2_ycoord,'FaceColor', [1 0 0]);
A6str = sprintf('Degree of cure = %6.3f', Deg_of_cure);
legend([a1],A6str)
end

spreadsheet_name = {source_files(i).name};
codes = {'Initial_peaktemp'    'Max_temperature'    'Final_peaktemp'    'Deg_of_cure'};
% maintable = array2table(zeros(0, numel(codes)+1), 'VariableNames', ['step', codes]);  %create empty destination table
% looptables = cell(len, 1);
looptable = table(spreadsheet_name,initial_peaktemp,Max_temperature,final_peaktemp,Deg_of_cure);
% looptable = table(spreadsheet_name,initial_peaktemp,Max_temperature,final_peaktemp,Deg_of_cure);
looptables{i} = looptable;
% filler = looptables{i},looptables{i+1};
maintable = vertcat(looptables{:});
writetable(maintable,'Results.xlsx');
end
%%-------------------------------------------------------------------------
%Find change in points
% pts=findchangepts(y,'Statistic','linear','MinThreshold',0.001);
%       plot(x,y,'-b',...
%           x(pts),y(pts),'rx')
