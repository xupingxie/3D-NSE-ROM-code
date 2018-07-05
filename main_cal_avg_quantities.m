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
function main_cal_avg_quantities(Draw_Models, quantity, fig_num, gen_elem_list, gen_avg_data)
global q_dim
tic
load(['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'],...
       'x','POD_u','POD_v','POD_w','center_u','center_v','center_w')
POD_u = POD_u'; center_u = center_u';
POD_v = POD_v'; center_v = center_v';
POD_w = POD_w'; center_w = center_w';
x_3d  = x;
% ---- on a z-plane
z_num_pts = 17;  z_pts = linspace(0.0625,  2.0625, z_num_pts); delta_z = z_pts(2)-z_pts(1);
x_num_pts = 100; x_pts = linspace(-7.4980, 7.4980, x_num_pts); delta_x = x_pts(2)-x_pts(1);
y_num_pts = 100; y_pts = linspace(-7.4980, 7.4980, y_num_pts); delta_y = y_pts(2)-y_pts(1);
% ---- choose model
Models      = {'DNS'; 'SMG'; 'VMS'; 'ML '; 'DS '; 'POD'};
if nargin==0
    Draw_Models = {'DNS'; ' '; ' '; 'ML '; ' '; 'POD'}; % determine what models are plotting
    quantity    = 'u';
    fig_num     = 1;
    gen_avg_data= 0;
end
cmp_res     = strcmp(Draw_Models, Models);
cmp_model   = []; for i=1:length(Models),...
        cmp_model = [cmp_model,num2str(cmp_res(i))];end
% ---- set time range
n_per_data = 4000;
n_per_start= 1;
n_per      = 4000;
T_set      = (n_per_start-1)*100+1:100:(n_per-1)*100+1;   %%%% change time interval %%%%
T_set_dns  = n_per_start:n_per;  T_set_smg  =T_set_dns;
delta_t    = 0.075;
skip       = 5;
str2       = []; tt = 0; legend_idx =[]; % record lengend info.

for i=1:length(Models) % ---- model
    if cmp_res(i)==1
        if i==1 % ---- for DNS
            str1 = '-*'; str2=[str2;'DNS        ']; tt=tt+1;
            load Matrices/r6/SNP_COEFF.dat
            d          = SNP_COEFF(:, 2:1+q_dim); a  = d';
            a          = a(:,T_set_dns);
            fprintf(1, 'DNS: \n');
        elseif i==2 % ---- for SMG
            str1 = 'r-o'; str2=[str2;'S-POD-ROM  ']; tt=tt+1;
            C_smg    = 500;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_hybrid_Rc4_Rz1_av'];
            Cur_name = [num2str(C_smg),'_',num2str(n_per_data), '_1_a'];
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set_smg);
            fprintf(1, 'SMG: \n');
        elseif i==3 % ---- for VMS
            str1 = 'g-x'; str2=[str2;'VMS-POD-ROM']; tt=tt+1;
            No_Large = 1;
            C_vms    = 1100; %Cav*10^6
            R_c      = 4;
            Alg      = 2;
            rare_upd = 1;
            Dir_load = ['Matrices/r', num2str(q_dim),'/connective_matrix145_193_17_VMS_'];
            Cur_name = ['L',num2str(No_Large),'hybrid_Rc', num2str(R_c), '_Rz1_av',num2str(C_vms),...
                '_', num2str(n_per_data),'_', num2str(rare_upd),'_alg',num2str(Alg),'_a']; % change this
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'VMS: \n');
        elseif i==4 % ---- for ML
            str1 = 'c-+'; str2=[str2;'ML-POD-ROM ']; tt=tt+1;
            C_ml     = 3000;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per_data),'_ML_C'];
            Cur_name = [num2str(C_ml),'_a'];
            Loadname = [Dir_load, Cur_name, '.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'ML : \n');
        elseif i==5 % ---- for DS
            str1 = 'm-s'; str2=[str2;'DS-POD-ROM ']; tt=tt+1;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145-193-17-',...
                num2str(n_per_data),'-DS-clip-n1.5e-01-a-every4-Rc4-newdelta'];
            Loadname = [Dir_load,'.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'DS : \n');
        elseif i==6 % ---- for POD-G
            str1 = 'k-d'; str2=[str2;'POD-G-ROM  ']; tt=tt+1;
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per_data),'_a'];
            Loadname = [Dir_load,'.mat'];
            load(Loadname, 'a');
            a        = a(:, T_set);
            fprintf(1, 'POD-G: \n');
        end
        
        savename   = ['Matrices/AVG/r',num2str(q_dim),'Avg_xzt_',quantity,'model_',num2str(i)];   % save data for avg
        
        if gen_avg_data == 1
            [N, M]     = size(a);
            snaps      = 1:skip:M;
            q_txyz           = zeros(length(snaps), x_num_pts, y_num_pts, z_num_pts);
            q_xavg           = zeros(length(snaps), y_num_pts, z_num_pts);
            q_xavg_zavg      = zeros(length(snaps), y_num_pts);
            q_xavg_zavg_tavg = zeros(y_num_pts,1);          
            % ---- interpolation
            for k = 1:z_num_pts
                z_k = z_pts(k);  % ---- on a z-plane
                index_zplane = find(x_3d(:,3)==z_k); % find pts on a z-plane
                x_2d         = x_3d(index_zplane,1:2); % get coordinates
                if mod(k,2)  == 1
                    Nx = 145; Ny = 193; ele_type = 'quadratic';
                else
                    Nx = 73;  Ny = 97;  ele_type = 'linear';
                end
                if (gen_elem_list==1) && (k==1||k==2)
                    first_time = 1; % generate mesh
                else
                    first_time = 0;
                end
                
                if first_time == 1
                    [x,e_conn,bnodes] = twod_mesh(0,1,0,1,'quadratic',Nx,Ny); % generate 2d mesh
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
                
                % fig_num = 1;
                % twod_plotm2(fig_num,x_2d,e_conn,'+')
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
                
                for t_m = 1:length(snaps);
                    mm   = snaps(t_m);
                    u    = center_u + POD_u*a(:,mm);
                    v    = center_v + POD_v*a(:,mm);
                    w    = center_w + POD_w*a(:,mm);
                    u_2d = u(index_zplane, 1);
                    v_2d = v(index_zplane, 1);
                    w_2d = w(index_zplane, 1);
                    
                    if strcmp(quantity, 'u')==1
                        q_2d = u_2d;
                    elseif strcmp(quantity, 'v')==1
                        q_2d = v_2d;
                    elseif strcmp(quantity, 'w')==1
                        q_2d = w_2d;
                    end
                    
                    for j = 1: y_num_pts % fix a y
                        y_j    = y_pts(j); % choose this number of points on the line y=y_j
                        for x_i = 1:x_num_pts
                        x_ydir = [x_pts(x_i) y_j];
                        q_txyz(t_m,x_i,j,k) = twod_interpolate(x_2d,e_conn,q_2d,x_ydir,...
                            element_list(x_i,j), bary1(x_i,j), bary2(x_i,j)); % u(t,x,y,z)
                        end
                    end                                        
                    Where_am_I( t_m,10)
                end
                Where_am_I(k,1)
            end
            % ---- take avg in x direction
            for k = 1:z_num_pts
                for t_m = 1:length(snaps);
                    for j = 1: y_num_pts % fix a y
                        q_xavg(t_m, j, k) = trapz(q_txyz(t_m,:,j,k))/(x_num_pts-1) ; % <u>_x
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
            
            save(savename, 'q_txyz', 'q_xavg', 'q_xavg_zavg', 'q_xavg_zavg_tavg');
        else
            load(savename)
        end
        
        figure(fig_num)
        legend_idx(1,tt)= plot(y_pts, q_xavg_zavg_tavg,str1,'linewidth',1.1);
        hold on
    end
    if i==length(Models)
        figure(fig_num)
        legend_handle= legend(legend_idx(1,:), str2);
        set(legend_handle,'FontSize',10)
        
        figure_avg = ['Matrices/AVG/figs/r',num2str(q_dim),'Avg_xzt_',quantity, cmp_model];
        xlabel('y'); ylabel(['<',quantity,'>_{x,z,t}'])
        set(gca,'FontSize',11)
        set(gca,'FontWeight','bold')
        print('-f', '-depsc',figure_avg)
    end
end
hold off
toc