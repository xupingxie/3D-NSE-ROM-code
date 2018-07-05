function [A,B,C] = offlane_compute_ROMmatrix(mu,filename,filename1,filename2,filenametu,filenametv,filenametw)

load(filename1)

POD_u = POD_u'; center_u = center_u';
POD_v = POD_v'; center_v = center_v';
POD_w = POD_w'; center_w = center_w';

[N         ,r      ] = size(POD_u);

load(filenametu)
for i=1:r
    %3-1
    cucux(i,1) = center_u'*Tu(i).ux*center_u;
    cupux(i,:) = center_u'*Tu(i).ux*POD_u;
    pucux(:,i) = POD_u'*Tu(i).ux*center_u;
    T(i).uu_x  = POD_u'*Tu(i).ux*POD_u;

    %3-2
    cvcuy(i,1) = center_v'*Tu(i).uy*center_u;
    cvpuy(i,:) = center_v'*Tu(i).uy*POD_u;
    pvcuy(:,i) = POD_v'*Tu(i).uy*center_u;
    T(i).vu_y  = POD_v'*Tu(i).uy*POD_u;

    %3-3
    cwcuz(i,1) = center_w'*Tu(i).uz*center_u;
    cwpuz(i,:) = center_w'*Tu(i).uz*POD_u;
    pwcuz(:,i) = POD_w'*Tu(i).uz*center_u;
    T(i).wu_z  = POD_w'*Tu(i).uz*POD_u;
end
clear Tu

load(filenametv)
for i=1:r
    %4-1
    cucvx(i,1) = center_u'*Tv(i).vx*center_v;
    cupvx(i,:) = center_u'*Tv(i).vx*POD_v;
    pucvx(:,i) = POD_u'*Tv(i).vx*center_v;
    T(i).uv_x  = POD_u'*Tv(i).vx*POD_v;

    %4-2
    cvcvy(i,1) = center_v'*Tv(i).vy*center_v;
    cvpvy(i,:) = center_v'*Tv(i).vy*POD_v;
    pvcvy(:,i) = POD_v'*Tv(i).vy*center_v;
    T(i).vv_y  = POD_v'*Tv(i).vy*POD_v;

    %4-3
    cwcvz(i,1) = center_w'*Tv(i).vz*center_v;
    cwpvz(i,:) = center_w'*Tv(i).vz*POD_v;
    pwcvz(:,i) = POD_w'*Tv(i).vz*center_v;
    T(i).wv_z  = POD_w'*Tv(i).vz*POD_v;
end
clear Tv

load(filenametw)
for i=1:r
    %5-1
    cucwx(i,1) = center_u'*Tw(i).wx*center_w;
    cupwx(i,:) = center_u'*Tw(i).wx*POD_w;
    pucwx(:,i) = POD_u'*Tw(i).wx*center_w;
    T(i).uw_x  = POD_u'*Tw(i).wx*POD_w;
    %5-2
    cvcwy(i,1) = center_v'*Tw(i).wy*center_w;
    cvpwy(i,:) = center_v'*Tw(i).wy*POD_w;
    pvcwy(:,i) = POD_v'*Tw(i).wy*center_w;
    T(i).vw_y  = POD_v'*Tw(i).wy*POD_w;
    %5-3
    cwcwz(i,1) = center_w'*Tw(i).wz*center_w;
    cwpwz(i,:) = center_w'*Tw(i).wz*POD_w;
    pwcwz(:,i) = POD_w'*Tw(i).wz*center_w;
    T(i).ww_z  = POD_w'*Tw(i).wz*POD_w;
end
clear Tw

for i=1:r
    C(i).mat = T(i).uu_x + T(i).vu_y + T(i).wu_z +...
        T(i).uv_x + T(i).vv_y + T(i).wv_z +...
        T(i).uw_x + T(i).vw_y + T(i).ww_z ;
end

load(filename2)
%3-4
cuxpux= center_u'*S11*POD_u;
Suxux = POD_u'*S11*POD_u;
%3-5
cuypuy= center_u'*S22*POD_u;
Suyuy = POD_u'*S22*POD_u;
%3-6
cvxpuy= center_v'*S12*POD_u;
Svxuy = POD_v'*S12*POD_u;
%3-7
cuzpuz= center_u'*S33*POD_u;
Suzuz = POD_u'*S33*POD_u;
%3-8
cwxpuz= center_w'*S13*POD_u;
Swxuz = POD_w'*S13*POD_u;

%4-4
cvxpvx= center_v'*S11*POD_v;
Svxvx = POD_v'*S11*POD_v;
%4-5
cuypvx= center_u'*S12'*POD_v;
Suyvx = POD_u'*S12'*POD_v;
%4-6
cvypvy= center_v'*S22*POD_v;
Svyvy = POD_v'*S22*POD_v;
%4-7
cvzpvz= center_v'*S33*POD_v;
Svzvz = POD_v'*S33*POD_v;
%4-8
cwypvz= center_w'*S23*POD_v;
Swyvz = POD_w'*S23*POD_v;

%5-4
cwxpwx= center_w'*S11*POD_w;
Swxwx = POD_w'*S11*POD_w;
%5-5
cuzpwx= center_u'*S13'*POD_w;
Suzwx = POD_u'*S13'*POD_w;
%5-6
cwypwy= center_w'*S22*POD_w;
Swywy = POD_w'*S22*POD_w;
%5-7
cvzpwy= center_v'*S23'*POD_w;
Svzwy = POD_v'*S23'*POD_w;
%5-8
cwzpwz= center_w'*S33*POD_w;
Swzwz = POD_w'*S33*POD_w;
clear S11 S12 S13 S22 S23 S33

%----zhu's original
A = cucux + cvcuy + cwcuz + 2*mu*(cuxpux') + mu*(cuypuy' + cvxpuy' + cuzpuz' + cwxpuz') +...
    cucvx + cvcvy + cwcvz + mu*(cvxpvx' + cuypvx') + 2*mu*(cvypvy') + mu*(cvzpvz' + cwypvz') +...
    cucwx + cvcwy + cwcwz + mu*(cwxpwx' + cuzpwx' + cwypwy' + cvzpwy') + 2*mu*cwzpwz' ;

B = cupux + pucux' + cvpuy + pvcuy' + cwpuz + pwcuz' + 2*mu*(Suxux') + mu*(Suyuy' + Svxuy' + Suzuz' + Swxuz') +...
    cupvx + pucvx' + cvpvy + pvcvy' + cwpvz + pwcvz' + mu*(Svxvx' + Suyvx') + 2*mu*(Svyvy') + mu*(Svzvz' + Swyvz') +...
    cupwx + pucwx' + cvpwy + pvcwy' + cwpwz + pwcwz' + mu*(Swxwx' + Suzwx' + Swywy' + Svzwy') + 2*mu*(Swzwz') ;

filename4 = [filename, 'GROM_ABC.mat'];
save(filename4, 'A','B','C')
