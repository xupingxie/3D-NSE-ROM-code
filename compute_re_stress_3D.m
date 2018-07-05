function compute_re_stress_3D(Draw_Models, quantity1, quantity2, fig_num, avg_q, gen_data, insetplot)
global q_dim
tic
q_dim = 6;
load(['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'],'x')
x_3d  = x;
% ---- on a z-plane
z_num_pts = 17 ; z_pts = linspace(0.0625,  2.0625, z_num_pts); delta_z = z_pts(2)-z_pts(1);
x_num_pts = 100; x_pts = linspace(-7.4980, 7.4980, x_num_pts); delta_x = x_pts(2)-x_pts(1);
y_num_pts = 100; y_pts = linspace(-7.4980, 7.4980, y_num_pts); delta_y = y_pts(2)-y_pts(1);
% ---- choose model
%Models      = {'DNS'; 'LDF'; 'EFROM'; 'PODG'; 'LProj'; 'ADROM' };
Models      = {'DNS'; 'LDF'; 'EFROM'; 'PODG'; 'LProj'; 'EFProj' };
if nargin==0
    Draw_Models = {'DNS'; 'LDF'; 'EFROM'; 'PODG'; 'LProj'; 'EFProj'}; % determine what models are plotting
    %Draw_Models = {'DNS'; 'LDF'; 'EFROM'; 'PODG'; 'LProj'; 'EFProj'}; % determine what models are plotting
%    Draw_Models = {'DNS'; ''; ''; 'PODG'; ''; 'ADROM'};
    quantity    = 'u';
    fig_num     = 1;
    avg_q       = 'tyz';
    quantity1='u';
    quantity2='v';
    insetplot=1;
    gen_data=0;
end

cmp_res     = strcmp(Draw_Models, Models);
cmp_model   = []; for i=1:length(Models),...
        cmp_model = [cmp_model,num2str(cmp_res(i))];end
% ---- set time range
n_per_data = 4000;
n_per_start= 1;
n_per      = 4000;
T_set      = (n_per_start-1)*100+1:100:(n_per-1)*100+1;   %%%% change time interval %%%%
T_set_dns  = n_per_start:n_per;  
delta_t    = 0.075;
skip       = 5;
str2       = []; tt = 0; legend_idx =[]; % record lengend info.

for i=1:length(Models) % ---- model
    %     if cmp_res(i)==1
    %         load(['Matrices/r',num2str(6),'/connective_matrix145_193_17.mat'],'x','POD_u','POD_v','POD_w',...
    %     'center_u','center_v','center_w')
    %     else
    %         load(['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'],'x','POD_u','POD_v','POD_w',...
    %     'center_u','center_v','center_w')
    %     end
    %     POD_u = POD_u'; center_u = center_u';
    %     POD_v = POD_v'; center_v = center_v';
    %     POD_w = POD_w'; center_w = center_w';
    snaps      = n_per_start:skip:n_per;
    if cmp_res(i)==1
        if i==1 % ---- for DNS
            str1 = '-'; str2=[str2;'DNS        ']; tt=tt+1;
            nameq='DNS';
            %             load Matrices/r6/SNP_COEFF.dat
            %             d          = SNP_COEFF(:, 2:1+q_dim); a  = d';
            %             a          = a(:,T_set_dns);
            fprintf(1, 'DNS: \n');
        elseif i==2 % ---- for LDF
            str1 = 'm-'; str2=[str2;'L-ROM-DF   ']; tt=tt+1;
            nameq='LDF-ROM';
            %             C_smg  = 500;
            %             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_hybrid_Rc4_Rz1_av'];
            %             Cur_name = [num2str(C_smg),'_',num2str(n_per_data), '_1_a'];
            %             Loadname = [Dir_load, Cur_name, '.mat'];
            %             load(Loadname, 'a');
            %             a        = a(:, T_set_smg);
            fprintf(1, 'LDF: \n');
        elseif i==5 % ---- for 
            str1 = 'g-'; str2=[str2;'L-ROM-Proj ']; tt=tt+1;
            nameq='LProj-ROM';
            %             No_Large = 1;
            %             C_vms    = 1100; %Cav*10^6
            %             R_c      = 4;
            %             Alg      = 2;
            %             rare_upd = 1;
            %             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_VMS_'];
            %             Cur_name = ['L',num2str(No_Large),'hybrid_Rc', num2str(R_c), '_Rz1_av',num2str(C_vms),...
            %                 '_', num2str(n_per_data),'_', num2str(rare_upd),'_alg',num2str(Alg),'_a']; % change this
            %             Loadname = [Dir_load, Cur_name, '.mat'];
            %             load(Loadname, 'a');
            %             a        = a(:, T_set);
            fprintf(1, 'LProj: \n');
        elseif i==3 % ---- for ML
            str1 = 'c-'; str2=[str2;'EF-ROM-DF  ']; tt=tt+1;
            nameq='EF-ROM';
            %             C_ml     = 3000;
            %             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
            %                 num2str(n_per_data),'_ML_C'];
            %             Cur_name = [num2str(C_ml),'_a'];
            %             Loadname = [Dir_load, Cur_name, '.mat'];
            %             load(Loadname, 'a');
            %             a        = a(:, T_set);
            fprintf(1, 'EF-ROM : \n');
        elseif i==6 % ---- for EF-Proj
            %str1 = 'r-'; str2=[str2;'AD-ROM     ']; tt=tt+1;
            %nameq='EFProj';
            %fprintf(1, 'EFProj : \n');
            
            str1 = 'r-'; str2=[str2;'EFProj-ROM ']; tt=tt+1;
            nameq   = 'EFProj';
            fprintf(1, 'EFProj : \n');
            
            
        elseif i==4 % ---- for POD-G
            str1 = 'k-'; str2=[str2;'G-ROM      ']; tt=tt+1;
            nameq='G-ROM';
            %             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
            %                 num2str(n_per_data),'_a'];
            %             Loadname = [Dir_load,'.mat'];
            %             load(Loadname, 'a');
            %             a        = a(:, T_set);
            fprintf(1, 'POD-G: \n');
        end
        
        
        if isempty(avg_q)~=1
            avgname   = ['Matrices/AVG/r',num2str(q_dim), 'Stress_',quantity1,quantity2,'_',avg_q,'_model_',nameq];
            %             % take average in different sense
            if gen_data == 1
                savename1   = ['Matrices/AVG/r',num2str(q_dim),'q_txyz_',nameq,'_',quantity1];    %----my stuff  
                savename2   = ['Matrices/AVG/r',num2str(q_dim),'q_txyz_',nameq,'_',quantity2];   
                %savename1   = ['Matrices/AVG/r',num2str(q_dim),'Avg_xzt_','umodel_',num2str(i)];    %--------Zhu's origin
                %savename2   = ['Matrices/AVG/r',num2str(q_dim),'Avg_xzt_','umodel_',num2str(i)];
                %savename1   = ['Matrices/AVG/r',num2str(q_dim),'Ptwise_',quantity1,'model_',num2str(i)];   % save data for avg
                %savename2   = ['Matrices/AVG/r',num2str(q_dim),'Ptwise_',quantity2,'model_',num2str(i)];   % save data for avg
                load(savename1)
                q1_txyz    = q_txyz;
                load(savename2)
                q2_txyz    = q_txyz;
                clear q_txyz
                
                avgname1    = ['Matrices/AVG/r',num2str(q_dim),'tyz_',nameq,'_',quantity1];    %--- my stuff    
                avgname2    = ['Matrices/AVG/r',num2str(q_dim),'tyz_',nameq,'_',quantity2];
                %avgname1    = ['Matrices/AVG/r',num2str(q_dim),'Avg_xzt_','umodel_',num2str(i)];            %---zhu's origin
                %avgname2    = ['Matrices/AVG/r',num2str(q_dim),'Avg_xzt_','umodel_',num2str(i)];
                %avgname1   = ['Matrices/AVG/r',num2str(q_dim),quantity1,'_',avg_q,'_model_',num2str(i)];
                %avgname2   = ['Matrices/AVG/r',num2str(q_dim),quantity2,'_',avg_q,'_model_',num2str(i)];
                load(avgname1)
                q1_avg_quantity = avg_quantity;
                load(avgname2)
                q2_avg_quantity = avg_quantity;        
                [points, stress]=stress_cal...
                    (q1_txyz,q2_txyz,q1_avg_quantity,q2_avg_quantity,snaps,skip,...
                    x_pts,y_pts,x_num_pts,y_num_pts,z_num_pts,avg_q);
                
                save(avgname, 'points', 'stress', 'quantity');
            else
                load(avgname)
            end
        end
        clear q1_txyz q2_txyz
        
        figure(fig_num)
        if insetplot==1
            
            if i~=4 && i~=6
                legend_idx(1,tt)= plot(points,stress,str1,'linewidth',1.1);
            else
                legend_idx(1,tt)= plot(points(1),stress(1),str1,'linewidth',1.1);
            end
            hold on
            
            if i==1||i==4||i==6     % plot DNS with EFProj and GROM
                figure(fig_num+5);
                plot(points, stress, str1,'linewidth',1.1);
                hold on
            end
            %xlabel('x','fontsize',15,'fontweight','b'); 
            %ylabel(['<  ',quantity1, '--< ',quantity1,' >, ', quantity2,'--< ', quantity2,' >  >'],'fontsize',14,'fontweight','b')
            set(gca,'FontSize',11)
            set(gca,'FontWeight','bold')
        else
            legend_idx(1,tt)= plot(points,stress,str1,'linewidth',1.1);
            hold on
        end
        
    end
    
    if isempty(avg_q)~=1
        str_t = ['_T',num2str(snaps(1)),'_',num2str(snaps(length(snaps)))];
        if i==length(Models)
            if strcmp(avg_q,'txz')
                x_axis = 'y';
            elseif strcmp(avg_q, 'tyz')
                x_axis = 'x';
            end
            figure(fig_num)
            legend_handle= legend(legend_idx(1,:), str2);
            set(legend_handle,'FontSize',10,'Position',[.75 .12 .15 .2])   %-----u-v
            %set(legend_handle,'FontSize',10,'Position',[.75 .7 .15 .2])   %-----u-w
            %set(legend_handle,'FontSize',10,'Position',[.6 .7 .15 .2])   %-----u-w
            figure_avg = ['Matrices/AVG/figs/r',num2str(q_dim), 'Stress_',avg_q,'_',quantity1,quantity2, cmp_model,str_t];
            %xlabel(x_axis,'fontsize',16,'fontweight','b'); 
            %ylabel(['<  ',quantity1, '--< ',quantity1,' >, ', quantity2,'--< ', quantity2,' >  >'],'fontsize',16,'fontweight','b')
            xlabel('$x$','interpreter','latex');
            ylabel('$\langle u - \langle u \rangle,v - \langle v \rangle \rangle$','interpreter','latex')
            set(gca,'FontSize',16)
            set(gca,'FontWeight','bold')
            %print('-f', '-depsc',figure_avg)
            
            % plot inset figures
             if insetplot ==1
                 fig1=figure(fig_num);
                 fig2=figure(fig_num+5);
                 %[h_m h_i]=inset(fig1,fig2);
                 inset(fig1,fig2);
             end
            figure(fig_num+5);
%            legend('DNS','G-ROM','AD-ROM')
%            legend_handle= legend(legend_idx(1,:), str2);
%            set(h,legend_handle,'FontSize',10,'Position',[.75 .12 .15 .2])
%            xlabel(x_axis,'fontsize',15,'fontweight','b'); 

            hold off
            
        end
    end
end
toc


function [points, stress]=stress_cal(q1,q2,q1_avg,q2_avg,snaps,skip,x_pts,y_pts,x_num_pts,y_num_pts,z_num_pts,avg_q)
temp = 1:skip:snaps(1)-1;  % to determin the starting position of t_m for snaps
if strcmp(avg_q,'txz')==1
    s_xavg           = zeros(length(snaps), y_num_pts, z_num_pts);
    s_xavg_zavg      = zeros(length(snaps), y_num_pts);
    s_xavg_zavg_tavg = zeros(y_num_pts,1);
    points = y_pts;
    
    % ---- calculate rms
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for j = 1: y_num_pts % fix a y
                s_xavg(t_m, j, k) = trapz(( q1(mm,:,j,k)-q1_avg(j,1) ).*( q2(mm,:,j,k)-q2_avg(j,1) ))...
                    /(x_num_pts-1) ; % <u>_x
            end
        end
    end
    % ---- take avg in z dirction
    for j = 1:y_num_pts
        for t_m = 1:length(snaps)
            s_xavg_zavg(t_m,j) = trapz(s_xavg(t_m,j,:))/(z_num_pts-1); %<u>_x,z
        end
    end
    % ---- take avg in time
    for j = 1:y_num_pts
        s_xavg_zavg_tavg(j,1) = trapz(s_xavg_zavg(:,j))/( length(snaps)-1 ); %<u>_x,z,t
    end
    stress = s_xavg_zavg_tavg;
    
elseif strcmp(avg_q,'tyz')==1
    s_yavg           = zeros(length(snaps), x_num_pts, z_num_pts);
    s_yavg_zavg      = zeros(length(snaps), x_num_pts);
    s_yavg_zavg_tavg = zeros(x_num_pts,1);
    
    points = x_pts;
    % ---- calculate rms
    % ---- take avg in y direction
    for k = 1:z_num_pts
        for t_m = 1:length(snaps);
            mm = t_m + length(temp);
            for x_i = 1: x_num_pts % fix a y
                s_yavg(t_m, x_i, k) = ...
                    trapz( ( q1(mm,x_i,:,k)-q1_avg(x_i,1) ).*( q2(mm,x_i,:,k)-q2_avg(x_i,1) )  )/(y_num_pts-1) ; % <u>_y
            end
        end
    end
    % ---- take avg in z dirction
    for x_i = 1:x_num_pts
        for t_m = 1:length(snaps)
            s_yavg_zavg(t_m,x_i) = trapz(s_yavg(t_m,x_i,:))/(z_num_pts-1); %<u>_y,z
        end
    end
    % ---- take avg in time
    for x_i = 1:x_num_pts
        s_yavg_zavg_tavg(x_i,1) = trapz(s_yavg_zavg(:,x_i))/( length(snaps)-1 ); %<u>_y,z,t
    end
    stress = s_yavg_zavg_tavg;
    
end
