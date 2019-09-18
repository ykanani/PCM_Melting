function [ Xi ] = zeroCrossing( x, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Find zero crossings:
Xi=nan(size(data));
cnv = data .* circshift(data,[1 0]);
slope = data - circshift(data,[1 0]);
[row, col] = find(cnv < 0 & abs(slope)>1e-5 );
% [row, col] = find(abs(slop(row1,col1)) > 1e-3 );
for i=1:size(data,2)
    zx=row(col==i);
    zx(zx==1)=[];
    zx(zx==size(data,1))=[];
    for k1 = 1:length(zx)
        i1=zx(k1)-1;
        i2=zx(k1);
        Xi(k1,i)= (0-data(i1,i))/(data(i2,i)-data(i1,i))*(x(i2)-x(i1))+x(i1);
        if Xi(k1,i)<0 
            aa=2;
        end
    end
    
end


% Interpolate zeros:

% figure
% plot(x, data)
% hold on
% % plot(an(zx), zeros(size(zx)), 'b*')
% plot(Xi, zeros(size(Xi)), 'pr')
% hold off
end

