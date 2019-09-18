function [Qtot Qlat]=calcQ(m,T,PCM)

Q=PCM.k*(T(end)-T(end-1))/m.Csize(end);
