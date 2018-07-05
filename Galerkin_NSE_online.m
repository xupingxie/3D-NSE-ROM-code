function dy = Galerkin_NSE_online(y,r,A,B,C)

dy = zeros(r,1);

for i=1:r
    dy(i) = -A(i,1) - B(i,:)*y - y'*C(i).mat*y;
end

return