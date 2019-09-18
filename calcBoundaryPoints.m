function Xfull=calcBoundaryPoints(X,method)
Xfull=X;
%%correct corner points
Xfull(1,1)=(Xfull(1,2)+Xfull(2,1))/2;
Xfull(end,1)=(Xfull(end,2)+Xfull(end-1,1))/2;

Xfull(1,end)=(Xfull(1,end-1)+Xfull(2,end))/2;
Xfull(end,end)=(Xfull(end,end-1)+Xfull(end-1,end))/2;

switch method
    case 'average'
        Xfull(1,:)=(Xfull(1,:)+Xfull(2,:))/2;
        Xfull(end,:)=(Xfull(end,:)+Xfull(end-1,:))/2;
        Xfull(:,1)=(Xfull(:,1)+Xfull(:,2))/2;
        Xfull(:,end)=(Xfull(:,end)+Xfull(:,end-1))/2;
    case 'zerograd'
        Xfull(1,:)=Xfull(2,:);
        Xfull(end,:)=Xfull(end-1,:);
        Xfull(:,1)=Xfull(:,2);
        Xfull(:,end)=Xfull(:,end-1);
end


