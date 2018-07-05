function adot = compute_derivative_coeff(a,dt)

[nrow, ncol]=size(a);
adot=zeros(nrow, ncol);
dh = 1/dt;

for j=1:ncol
    if j==1
       tmp=(-25/12)*a(:,j)+4*a(:,j+1)-3*a(:,j+2)+(4/3)*a(:,j+3)-(1/4)*a(:,j+4);
       adot(:,j)= tmp*dh;
    elseif j==2
       tmp=(-25/12)*a(:,j)+4*a(:,j+1)-3*a(:,j+2)+(4/3)*a(:,j+3)-(1/4)*a(:,j+4);
       adot(:,j)= tmp*dh;
    elseif j==ncol-1
       tmp=(25/12)*a(:,j)-4*a(:,j-1)+3*a(:,j-2)-(4/3)*a(:,j-3)+(1/4)*a(:,j-4);
       adot(:,j)= tmp*dh;
    elseif j==ncol
       tmp=(25/12)*a(:,j)-4*a(:,j-1)+3*a(:,j-2)-(4/3)*a(:,j-3)+(1/4)*a(:,j-4);
       adot(:,j)= tmp*dh;
    else
       tmp=(1/12)*a(:,j-2)-(2/3)*a(:,j-1)+(2/3)*a(:,j+1)-(1/12)*a(:,j+2);
       adot(:,j)= tmp*dh;
 end
 
 end
    