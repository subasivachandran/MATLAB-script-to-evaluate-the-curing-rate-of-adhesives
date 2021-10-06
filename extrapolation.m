% curve = xlsread('curve.xlsx');
% day = curve(:,1);
% hairBalls = curve(:,2);
% steps2 = (-5:0.2:25);
% extrapoSpline2 = interp1(day,hairBalls,steps2,'spline','extrap');
% figure(2)
% plot(day, hairBalls,'ko', steps2,extrapoSpline2,'gx:')
% % steps = (0:.2:20);
% % interSpline = interp1(day,hairBalls,steps,'spline');
% % extrapoSpline = interp1(day,hairBalls,steps,'spline','extrap');
% % plot(day, hairBalls,'ko',steps,interSpline,'rx:');

% A=[1 2 3 2 1 4 5 6 7 5 2 0 4 5 6 10]
% [peaks,idx]=findpeaks(A)


% A=[5 2 3 2 1 4 5 6 7 5 2 0 4 5 6 10]
% d=sign([0 diff(A)]);
% idx=strfind(d,[1,-1])

% A = [1 2 1 2 1 2 3 4 5 6 7 5 8 0 1 3 6 9 4];
% da = diff(A)>0
% stats = regionprops(bwlabel(da), 'Area', 'PixelIdxList')
% % Find which runs have a length of at least
% % 3 elements of monotonically increasing elements
% threeOrLonger = find([stats.Area] >= 3)
% % Print out those sequences.
% if ~isempty(threeOrLonger)
%   for blobIndex = 1 : length(threeOrLonger)
%     % This blob is 3. Find where the blob starts and stops.
%     startingIndex = stats(threeOrLonger(blobIndex)).PixelIdxList(1);
%     endingIndex = stats(threeOrLonger(blobIndex)).PixelIdxList(end)+1;
%     % Print out that run.
%     fprintf('For run #%d, startingIndex = %d, endingIndex = %d, the elements = ', ...
%       blobIndex, startingIndex, endingIndex);
%     fprintf('%d', A(startingIndex : endingIndex));
%     fprintf('\n');
%   end
% else
%   msgbox('No run is 3 or longer');
% end
% data = readtable('testdata.xlsx');
% x = data.Temp;
% y = data.DSC;
% pts=findchangepts(y,'Statistic','linear','MinThreshold',0.01);
%       plot(x,y,'-b',...
%           x(pts),y(pts),'rx')
% y=[1 2 5 0 4 2 7 8 9]
%  max([1 max((find(diff(y)<=0))+1)])

% function [y1, y2, y3] = MyInterpolator()
% x =[ 1 2  3 4  5  6   7  8  9  10]
% y=[ -2 3 -4 6 -7 10 -17 25 -26 30]
% plot(x,y,'ob')
% xi=1:0.1:10
% method=char('linear','cubic','spline')
% col={'--r',':g','-k'}
% for k=1:3
% yi{k}=interp1(x,y,xi,method(k,:))
% hold on;plot(xi,yi{k},col{k})
% end
% y1=yi{1}
% y2=yi{2}
% y3=yi{3}
% 
% x = linspace(0, 50);                                            % Create Data
% y = 6 * (exp(-(x-10).^2/25) + exp(-(x-25).^2/50));              % Create Data
% Lvy = (y >  1) & (x > 20);
% figure
% plot(x, y)
% hold on
% patch([x(Lvy) fliplr(x(Lvy))], [ones(size(x(Lvy))) fliplr(y(Lvy))], 'g')
% hold off

data = readtable('testdata.xlsx');
x = data.Temp;
y = data.DSC; 
y_0 = y*0;
% x=0:0.1:3.14;
% y=sin(4*x).*sec(3*x);
xc = find((circshift(y, [0 1])) > 0); % Define ?y? Zero-Crossings
length_xc = length(xc);
idx_begin = xc(1)-1:xc(1);% Index Range Including Zero-Crossing
idx_end = xc(length_xc):xc(length_xc)+1;
x_begin = interp1(y(xc(1)-1:xc(1)), x(xc(1)-1:xc(1)), 0);        % Find ?x? At Zero Crossing
x_end = interp1(y(xc(length_xc):xc(length_xc)+1), x(xc(length_xc):xc(length_xc)+1), 0);
plot(x,y,x_begin,0,'x',x_end,0,'x',x,y_0,'-');