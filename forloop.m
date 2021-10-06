% clear all;
% clc;
% x = [1:1:10];
% y = [11:1:20];
% m = numel(x);
% variables_new = [];
% [gy] = gradient(y);
% [gx] = gradient(x);
% for k=1:m
%    % insert your if-else code here to test x(k) for some condition and set g(k)
%    if x(k) >= 5 & x(k) <= 10
%        variables_temp = [(gy(k) ./ gx(k)),x(k)]
%        variables_new = [variables_new;variables_temp]
%    end
% end

% x= [90.3384900000000,91.8425100000000,97.2637500000000,135.963910000000,90.3384900000000];
% y = [-0.0964600000000000,0,0,-0.100390000000000,-0.0964600000000000];
% plot(x,y)
% fill(x,y,'r')
% hold on
% % a3 = area(x,y,'FaceColor','g');
% % A3str = sprintf('Area 3 = %6.3f', a3)
% a4 = polyarea(x,y);
% a4_trapz = trapz(x,y);
% afinal = a4 - A3str;
% % afinal = area(x,y,'FaceColor','g');
% Afinal = sprintf('Area final = %6.3f', a4_trapz);
% legend([a3],A3str)
% hold off

% x= [90.616840000000000,93.813606121842510,98.837425198610530,1.538424800000000e+02,90.616840000000000]
% y= [0.268280000000000,0,0,0.263970000000000,0.268280000000000]
% poar = trapz(x,y);
% plot(x,y)

%Linspace
% startvalue = -1.3750;
% endvalue = -1.2013;
% numberofpoints = 688; %includes start and end point
% values = linspace(startvalue, endvalue, numberofpoints);

% x=[0:0.1:10];
% x1=5;
% plot(x,2*sin(x),x,5*sin(x),x,cos(x));
% y1=get(gca,'ylim')
% hold on
% plot([x1 x1],y1)
codes = {'ABO1'    'ABO2'    'ADE3'    'ADE4'    'AGUG' 'ACE4'};
maintable = array2table(zeros(0, numel(codes)+1), 'VariableNames', ['step', codes]);  %create empty destination table
%your loop replace by actual code
% numsteps = 10;
% looptables = cell(numsteps, 1);
% for step = 1:numsteps   %your loop, replace by actual code
%   %random generation of loop table, replace by actual code
%   randrows = randi(numel(codes));
%   looptable = table(codes(randperm(numel(codes), randrows))', rand(randrows, 1)*20000, 'VariableNames', {'codes', 'values'});
%   looptables{step} = looptable;
% %     %filling of the main table
% %     filler = nan(1, numel(codes));
% %     [~, destcol] = ismember(looptable.codes, codes);
% %     filler(destcol) = looptable.values;
% %     maintable{step, :} = [step, filler];
% maintable2 = vertcat(looptables{:});
% maintable2.step = repelem((1:numel(looptables))', cellfun(@height, looptables));
% maintable2 = unstack(maintable2, 'values', 'codes')
%   end

T = table([10;20;30],{'M';'F';'F'},'VariableNames',{'Age','Gender'},'RowNames',{'P1','P2','P3'})
