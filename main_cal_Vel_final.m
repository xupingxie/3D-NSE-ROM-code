% Zhu Wang 07-21-2011
% Codes for generating new figures for our POD-comparison paper.
%
% to calculate < >_x or < >_y, we need finite element interpolation,
% here, we employ Dr. Borggaard's codes, do 2d-interpolation on each
% z-plane.
% Since we skip the center of the cube when we generate 3d mesh (one
% cube is split into five tetrahedrons), in z-direction, on z-plane,
% we re-construct triangular meshes for quadratic, linear, quadratic,
% linear, ... elements one by one.
%
% domain x: [-7.4980, 7.4980]; 145 pts on a radius
%        y: [-7.4980, 7.4980]; 192 pts in azimuthal direction
%        z: [0.0625,  2.0625]; 17  pts in spanwise direction
% z=0.0625, num of nodes 192*145=27840, recover mesh for quadratic elements
% z=0.1875, num of nodes 96*145+96*73=20928, recover mesh for linear elements
% z=0.3125, num of nodes 27840, ...
% ...
%
% For example, <u>_{x,z,t} is a function of y.
% step 1, fix Z_i and construct triangular mesh (quadratic or linear);
%         set u
% step 2, fix y_j, find u(x_k,y_j) for k=1:N by 2d interpolation;
% step 3, do time average
% step 4, plot
%
% generate u to be a 4-dim array; t-h, x-i, y-j, z-k
function main_cal_Vel_final(Draw_Models, quantity, fig_num, gen_elem_list, gen_data, avg_q, load_avg)
global q_dim
q_dim = 6;
tic
load(['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'],'x')
% POD_u = POD_u'; center_u = center_u';
% POD_v = POD_v'; center_v = center_v';
% POD_w = POD_w'; center_w = center_w';
x_3d  = x;
% ---- on a z-plane
z_num_pts = 17 ; z_pts = linspace(0.0625,  2.0625, z_num_pts); delta_z = z_pts(2)-z_pts(1);
x_num_pts = 100; x_pts = linspace(-7.4980, 7.4980, x_num_pts); delta_x = x_pts(2)-x_pts(1);
y_num_pts = 100; y_pts = linspace(-7.4980, 7.4980, y_num_pts); delta_y = y_pts(2)-y_pts(1);
% ---- choose model
Models      = {'DNS'; 'POD'; 'EFROM'; 'EFProj'; 'LDF'; 'LProj' };
if nargin==0
    Draw_Models = {'DNS'; 'POD'; 'EFROM'; 'EFProj'; 'LDF'; 'LProj' }; % determine what models are plotting
    quantity    = 'w';
    fig_num     = 1;
    gen_elem_list = 0;
    gen_data    = 1;
    avg_q       = 'tyz';  %[]
    load_avg    = 1;
    insetplot   = 1;
end
cmp_res     = strcmp(Draw_Models, Models);
cmp_model   = []; for i=1:length(Models),...
        cmp_model = [cmp_model,num2str(cmp_res(i))];end
% ---- set time range
n_per_data = 4000;
n_per_start= 1;
n_per      = 4000;
T_set      = (n_per_start-1)*100+1:100:(n_per-1)*100+1;   %%%% change time interval %%%%
T_set_dns  = n_per_start:n_per;  T_set_smg  =T_set;%T_set_DNS
delta_t    = 0.075;
skip       = 5;
str2       = []; tt = 0; legend_idx =[]; % record lengend info.

for i=1:length(Models) % ---- model
    if cmp_res(1)==1
        load(['Matrices/r',num2str(6),'/connective_matrix145_193_17.mat'],'x','POD_u','POD_v','POD_w',...
            'center_u','center_v','center_w')
    else
        load(['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'],'x','POD_u','POD_v','POD_w',...
            'center_u','center_v','center_w')
    end
    POD_u = POD_u'; center_u = center_u';
    POD_v = POD_v'; center_v = center_v';
    POD_w = POD_w'; center_w = center_w';
     
    if cmp_res(i)==1
        if i==1 % ---- for DNS
            str1 = '-'; str2=[str2;'DNS        ']; tt=tt+1;
            load Matrices/r6/SNP_COEFF.dat
            nameq = 'DNS';
            d          = SNP_COEFF(:, 2:7); a  = d';
            a          = a(:,T_set_dns);
            fprintf(1, 'DNS: \n');
        elseif i==4 % ---- for EFProj
            str1 = 'r-'; str2=[str2;'EF-ROM-Proj']; tt=tt+1;
            r1 = 4;
            nameq   = 'EFProj';
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_'];
            Cur_name = [num2str(n_per_data), 'EFProj_',num2str(r1),'_a'];
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'EFProj: \n');
% %             str1 = 'r-'; str2=[str2;'AD-ROM     ']; tt=tt+1;
% %             %r1 = 4;
% %             de=0.3;
% %             nu=0.0285;
% %             nameq   = 'ADROM';
% %             %Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_'];
% %             %Cur_name = [num2str(n_per_data), 'EFProj_',num2str(r1),'_a'];
% %             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
% %                 num2str(n_per),'ADROM_new_FEDF'];
% %             Loadname = [Dir_load, num2str(de),'_',num2str(nu),'_a.mat'];    
% %             %Cur_name = ['ADROM_new_FEDF'] ;
% %             %Loadname = [Dir_load, Cur_name, '.mat'];
% %             load(Loadname, 'a');
% %             a        = a(:, T_set);
% %            fprintf(1, ': \n');
        elseif i==5 % ---- for LDF
            str1 = 'm-'; str2=[str2;'L-ROM-DF   ']; tt=tt+1;
            nameq='LDF-ROM';
            de = 0.247;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_'];
            Cur_name = [num2str(n_per_data), 'LDF_', num2str(de),'_a'];
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'LDF: \n');
        elseif i==3 % ---- for EFROM
            str1 = 'c-'; str2=[str2;'EF-ROM-DF  ']; tt=tt+1;
            nameq='EF-ROM';
            de        = 0.001367;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per_data),'EFROM_new'];
            Cur_name = [num2str(de),'_a'];
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'EFROM : \n');
        elseif i==6 % ---- for LProj
            str1 = 'g-'; str2=[str2;'L-ROM-Proj ']; tt=tt+1;
            r1        = 1;
            nameq='LProj-ROM';
            Dir_load = ['Matrices/r', num2str(q_dim),'/connective_matrix145_193_17_'];
            Cur_name = [num2str(n_per_data), 'LProj_', num2str(r1), '_a']; % change this
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'LProj: \n');
        elseif i==2 % ---- for POD-G
            str1 = 'k-'; str2=[str2;'G-ROM      ']; tt=tt+1;
            nameq='G-ROM'
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per_data),'_GROM_a'];
            Loadname = [Dir_load,'.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'POD-G: \n');
        end
        
        %savename   = ['Matrices/AVG/r',num2str(q_dim),'Ptwise_',quantity,'model_',num2str(i)];   % save data for avg
        savename   = ['Matrices/AVG/r', num2str(q_dim), 'q_txyz_', nameq, '_', quantity];   % save data for avg
        [N, M]     = size(a);
        snaps_full = 1:skip:M;
        snaps      = n_per_start:skip:n_per;
        
        if load_avg ~=1
            if gen_data == 1
                q_txyz           = zeros(length(snaps_full), x_num_pts, y_num_pts, z_num_pts);
                % ---- interpolation
                for k = 1:z_num_pts
                    z_k = z_pts(k);  % ---- on a z-plane
                    index_zplane = find(x_3d(:,3)==z_k); % find pts on a z-plane
                    if mod(k,2)  == 1
                        Nx = 145; Ny = 193; ele_type = 'quadratic';
                    else
                        Nx = 73;  Ny = 97;  ele_type = 'linear';
                        for j=1:96
                            index_skip((j-1)*73+1:j*73) = (j-1)*(73+145)+(1:2:145);
                        end
                        index_skip(96*73+1:97*73) = index_skip(1:73);
                        index_zplane = index_zplane(index_skip);
                    end
                    x_2d         = x_3d(index_zplane,1:2); % get coordinates
                    
                    if (gen_elem_list==1) && (k==1||k==2)
                        first_time = 1; % generate mesh
                    else
                        first_time = 0;
                    end
                    
                    if first_time == 1
                        [x,e_conn,bnodes] = twod_mesh(0,1,0,1,ele_type,Nx,Ny); % generate 2d mesh
                        clear x bnodes % clear coordinates and boundary, just save connective matrix
                        index        = 1:Nx*Ny;
                        index( Nx*(Ny-1)+1 : Nx*Ny ) = 1:Nx; % two radii coincide with each other
                        e_conn       = index( e_conn ); % change connective matrix correspondingly, 145*192+1-> 1
                        %     neighbors = tri_mesh_neighbors_new( e_conn );
                        [n_elem,~]   = size(e_conn);
                        barycenter = zeros(n_elem, 2);
                        for x_i= 1:n_elem
                            barycenter(x_i,:) = [sum(x_2d(e_conn(x_i,1:3),1))/3, sum(x_2d(e_conn(x_i,1:3),2))/3];
                        end
                        save(['Matrices/r6/Neighbor_info_',num2str(Nx),'_',num2str(Ny')],'e_conn','barycenter')
                    else
                        load(['Matrices/r6/Neighbor_info_',num2str(Nx),'_',num2str(Ny')]);
                    end
                    
                    %                 fig_num = 1;
                    %                 twod_plotm2(fig_num,x_2d,e_conn,'+')
                    % ---- generate file that list the element that a point locates
                    elem_name = ['Matrices/AVG/',ele_type,'_elementlist.mat'];
                    if first_time == 1
                        element_list = zeros(x_num_pts, y_num_pts);
                        bary1        = zeros(x_num_pts, y_num_pts);
                        bary2        = zeros(x_num_pts, y_num_pts);
                        for j = 1: y_num_pts % fix a y
                            y_j    = y_pts(j); % choose this number of points on the line y=y_j
                            x_ydir = [x_pts' y_j*ones(x_num_pts,1)];
                            for x_i = 1: x_num_pts
                                xy_point = x_ydir(x_i,:);
                                [element_list(x_i,j), bary1(x_i,j), bary2(x_i,j)] ...
                                    = mesh_search(barycenter, x_2d, e_conn(:,1:3), xy_point);
                            end
                        end
                        save(elem_name,'element_list','bary1','bary2');
                    else
                        load(elem_name);
                    end
                    
                    disp('compute q_txyz now')
                    for t_m = 1:length(snaps_full);
                        mm   = snaps_full(t_m);
                        
                        if strcmp(quantity, 'u')==1
                            u    = center_u + POD_u*a(:,mm);
                            u_2d = u(index_zplane, 1);
                            q_2d = u_2d;
                        elseif strcmp(quantity, 'v')==1
                            v    = center_v + POD_v*a(:,mm);
                            v_2d = v(index_zplane, 1);
                            q_2d = v_2d;
                        elseif strcmp(quantity, 'w')==1
                            w    = center_w + POD_w*a(:,mm);
                            w_2d = w(index_zplane, 1);
                            q_2d = w_2d;
                        end
                        
                        for j = 1: y_num_pts % fix a y
                            y_j     = y_pts(j); % choose this number of points on the line y=y_j
                            for x_i = 1:x_num_pts
                                x_ydir  = [x_pts(x_i) y_j];
                                q_txyz(t_m,x_i,j,k) = twod_interpolate(x_2d,e_conn,q_2d,x_ydir,...
                                    element_list(x_i,j), bary1(x_i,j), bary2(x_i,j)); % u(t,x,y,z)
                            end
                        end
                        %Where_am_I( t_m,10)
                        sprintf('t_m= %d,k= %d ', t_m, k)
                    end
                    %Where_am_I(k,1)
                end
                save(savename, 'q_txyz');
            else
                load(savename)
            end
        end
        if isempty(avg_q)~=1
            %avgname   = ['Matrices/AVG/r',num2str(q_dim),quantity,'_',avg_q,'_model_',num2str(i)];
            avgname   = ['Matrices/AVG/r',num2str(q_dim),avg_q,'_',nameq,'_',quantity];
            
            % take average in different sense
            if load_avg == 1
                load(avgname);
            else
                [points, avg_quantity, rms_quantity]=take_avg...
                    (q_txyz,avg_q,snaps,skip,x_pts,y_pts,x_num_pts,y_num_pts,z_num_pts);
                save(avgname, 'points', 'avg_quantity', 'rms_quantity', 'quantity');
            end
            
            %figure(fig_num)
            %legend_idx(1,tt)= plot(points, avg_quantity,str1,'linewidth',1.1);
            %hold on
            figure(fig_num+20)
            legend_idx2(1,tt)= plot(points, rms_quantity,str1,'linewidth',1.1);
            hold on
            
            figure(fig_num)
            if insetplot==1
                
                if i~=2 && i~=4
                    legend_idx(1,tt)= plot(points,avg_quantity,str1,'linewidth',1.1);
                else
                    legend_idx(1,tt)= plot(points(1),avg_quantity(1),str1,'linewidth',1.1);
                end
                hold on
            
                %if i==1||i==3||i==6||i==5
                %    legend_idx(1,tt)= plot(points, avg_quantity,str1,'linewidth',1.1);
                %    hold on
                %end
                
                figure(fig_num+5)
                if i==1||i==2||i==4     % plot DNS with EFProj and GROM
                    figure(fig_num+5);
                    plot(points, avg_quantity,str1,'linewidth',1.1);
                    hold on
                end
                %set(gca,'FontSize',11)
                %set(gca,'FontWeight','bold')
%                 legend_idx3(1,tt)= plot(points, avg_quantity,str1,'linewidth',1.1);
%                 hold on
                %legend_idx(1,tt)= plot(points, avg_quantity,str1,'linewidth',1.1);
      
%                 if i~=4 && i~=6
%                     legend_idx(1,tt)= plot(points,stress,str1,'linewidth',1.1);
%                 else
%                     legend_idx(1,tt)= plot(points(1),stress(1),str1,'linewidth',1.1);
%                 end
%                 hold on
%                 
%                 if i==1||i==4||i==6     % plot DNS with EFProj and GROM
%                     figure(fig_num+5);
%                     plot(points, stress, str1,'linewidth',1.1);
%                     hold on
%                 end
                %xlabel('x','fontsize',15,'fontweight','b');
                %ylabel(['<  ',quantity1, '--< ',quantity1,' >, ', quantity2,'--< ', quantity2,' >  >'],'fontsize',14,'fontweight','b')
%                 set(gca,'FontSize',11)
%                 set(gca,'FontWeight','bold')
            end
        
        end
        
    end
    if i==length(Models)&&sum(cmp_res)>0
        str_t = ['_T',num2str(snaps(1)),'_',num2str(snaps(length(snaps)))];
        if strcmp(avg_q,'txz')
            x_axis = 'y';
        elseif strcmp(avg_q, 'tyz')
            x_axis = 'x';
        end
        figure(fig_num)
        legend_handle= legend(legend_idx(1,:), str2);
        set(legend_handle,'FontSize',10)
        
        figure_avg = ['Matrices/AVG/figs/r',num2str(q_dim),...
            'Avg_',avg_q,'_',quantity, cmp_model,str_t];
        %xlabel(x_axis); ylabel(['<',quantity,'>']);
        xlabel('$x$','interpreter','latex');
        ylabel('$\langle w \rangle$','interpreter','latex')
        set(gca,'FontSize',16)
        set(gca,'FontWeight','bold')
        %print('-f', '-depsc',figure_av
        hold off
        
        %figure(fig_num*10)
        if insetplot ==1
            fig1=figure(fig_num);
            fig2=figure(fig_num+5);
            %[h_m h_i]=inset(fig1,fig2);
            inset(fig1,fig2);
        end
     
        
        
        figure(fig_num+20)
        legend_handle= legend(legend_idx2(1,:), str2);
        set(legend_handle,'FontSize',10)
        
        figure_avg = ['Matrices/AVG/figs/r',num2str(q_dim),...
            'rms_',avg_q,'_',quantity, cmp_model,str_t];
        %xlabel(x_axis); ylabel(['<',quantity,'>_{rms}'])
        xlabel('$x$','interpreter','latex');
        ylabel('$\langle w \rangle_{rms}$','interpreter','latex');
        set(gca,'FontSize',16)
        set(gca,'FontWeight','bold')
        %print('-f', '-depsc',figure_avg)
       
        hold off
        
    end
end
%
toc

function [points, avg_quantity, rms_quantity]=take_avg(q_txyz,avg_q,snaps,skip,x_pts,y_pts,x_num_pts,y_num_pts,z_num_pts)
temp = 1:skip:snaps(1)-1;  % to determin the starting position of t_m for snaps
if strcmp(avg_q,'txz')==1
    q_xavg           = zeros(length(snaps), y_num_pts, z_num_pts);
    q_xavg_zavg      = zeros(length(snaps), y_num_pts);
    q_xavg_zavg_tavg = zeros(y_num_pts,1);
    q_xrms           = zeros(length(snaps), y_num_pts, z_num_pts);
    q_xrms_zrms      = zeros(length(snaps), y_num_pts);
    q_xrms_zrms_trms = zeros(y_num_pts,1);
    % ---- take avg in x direction
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for j = 1: y_num_pts % fix a y
                q_xavg(t_m, j, k) = trapz(q_txyz(mm,:,j,k))/(x_num_pts-1) ; % <u>_x
            end
        end
    end
    % ---- take avg in z dirction
    for j = 1:y_num_pts
        for t_m = 1:length(snaps)
            q_xavg_zavg(t_m,j) = trapz(q_xavg(t_m,j,:))/(z_num_pts-1); %<u>_x,z
        end
    end
    % ---- take avg in time
    for j = 1:y_num_pts
        q_xavg_zavg_tavg(j,1) = trapz(q_xavg_zavg(:,j))/( length(snaps)-1 ); %<u>_x,z,t
    end
    points = y_pts; avg_quantity = q_xavg_zavg_tavg;
    
    % ---- calculate rms
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for j = 1: y_num_pts % fix a y
                q_xrms(t_m, j, k) = trapz(( q_txyz(mm,:,j,k)-avg_quantity(j,1) ).^2)...
                    /(x_num_pts-1) ; % <u>_x
            end
        end
    end
    % ---- take avg in z dirction
    for j = 1:y_num_pts
        for t_m = 1:length(snaps)
            q_xrms_zrms(t_m,j) = trapz(q_xrms(t_m,j,:))/(z_num_pts-1); %<u>_x,z
        end
    end
    % ---- take avg in time
    for j = 1:y_num_pts
        q_xrms_zrms_trms(j,1) = trapz(q_xrms_zrms(:,j))/( length(snaps)-1 ); %<u>_x,z,t
    end
    rms_quantity = sqrt(q_xrms_zrms_trms);
    
elseif strcmp(avg_q,'tyz')==1
    q_yavg           = zeros(length(snaps), x_num_pts, z_num_pts);
    q_yavg_zavg      = zeros(length(snaps), x_num_pts);
    q_yavg_zavg_tavg = zeros(x_num_pts,1);
    q_yrms           = zeros(length(snaps), x_num_pts, z_num_pts);
    q_yrms_zrms      = zeros(length(snaps), x_num_pts);
    q_yrms_zrms_trms = zeros(x_num_pts,1);
    % ---- take avg in y direction
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for x_i = 1: x_num_pts % fix a y
                q_yavg(t_m, x_i, k) = ...
                    trapz(q_txyz(mm,x_i,:,k))/(y_num_pts-1) ; % <u>_y
            end
        end
    end
    % ---- take avg in z dirction
    for x_i = 1:x_num_pts
        for t_m = 1:length(snaps)
            q_yavg_zavg(t_m,x_i) = trapz(q_yavg(t_m,x_i,:))/(z_num_pts-1); %<u>_y,z
        end
    end
    % ---- take avg in time
    for x_i = 1:x_num_pts
        q_yavg_zavg_tavg(x_i,1) = trapz(q_yavg_zavg(:,x_i))/( length(snaps)-1 ); %<u>_y,z,t
    end
    points = x_pts; avg_quantity = q_yavg_zavg_tavg;
    
    % ---- calculate rms
    % ---- take avg in y direction
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for x_i = 1: x_num_pts % fix a y
                q_yrms(t_m, x_i, k) = ...
                    trapz(( q_txyz(mm,x_i,:,k)-avg_quantity(x_i,1) ).^2  )/(y_num_pts-1) ; % <u>_y
            end
        end
    end
    % ---- take avg in z dirction
    for x_i = 1:x_num_pts
        for t_m = 1:length(snaps)
            q_yrms_zrms(t_m,x_i) = trapz(q_yrms(t_m,x_i,:))/(z_num_pts-1); %<u>_y,z
        end
    end
    % ---- take avg in time
    for x_i = 1:x_num_pts
        q_yrms_zrms_trms(x_i,1) = trapz(q_yrms_zrms(:,x_i))/( length(snaps)-1 ); %<u>_y,z,t
    end
    rms_quantity = sqrt(q_yrms_zrms_trms);
    
end