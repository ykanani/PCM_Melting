function Xavg=volAvg(m,X)
Xavg=sum(sum(X(2:m.Nx+1,2:m.Ny+1).*m.V(2:m.Nx+1,2:m.Ny+1)))/sum(sum(m.V(2:m.Nx+1,2:m.Ny+1)));

end
