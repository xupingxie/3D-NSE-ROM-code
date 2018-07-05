% This subroutin is for plotting Energy Spectrum;
% Call subroutine FFT_EgySpm.m
clear
global q_dim
q_dim  = 6;
%Models = {'DNS','ADL','EFProj','LDF','LProj','GROM'};
Models = {'DNS','ADL','ED','LDF','LProj','GROM'};
Draw_Models = {'DNS','','','','','GROM'}; % determine what models are plotting
Load_FFT =[1,0,0,0,0,0]; % 1- load FFT results for corresponding model

cmp_res = strcmp(Draw_Models, Models);
cmp_model = []; for i=1:length(Models),...
        cmp_model = [cmp_model,num2str(cmp_res(i))];end
delta_t = 0.075;

str2 =[];
tt   = 0; legend_idx =[]; % record lengend info.
do_fig =1; % not plot energy spectrum for each model
POD_data = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17.mat'];
load(POD_data)
%---------
do_avg   = 1; % choose do_avg to be 1 and make spatial
% avg to draw egy spectrum.
x_cord = [0.9992, 0.3575, 1.0625]; % Probe/ center of Probes
if do_avg==0
    Probe_index = find(abs(x(:,1)-x_cord(1))<0.001&abs(x(:,2)-x_cord(2))<0.001...
        &abs(x(:,3)-x_cord(3))<0.001); % just one point
    str_avg = '_1P';
else
    Probe_index = find(abs(x(:,1)-x_cord(1))<0.1...
        &abs(x(:,2)-x_cord(2))<0.1...
        &abs(x(:,3)-x_cord(3))<0.2); % choose a box then take avg.
    str_avg = ['_',num2str(length(Probe_index)), 'P'];
end
P_cord = x(Probe_index,:);
fprintf(1, [num2str(min(P_cord(:,1))),' , ', num2str(max(P_cord(:,1))), ' \n']);
fprintf(1, [num2str(min(P_cord(:,2))),' , ',  num2str(max(P_cord(:,2))), ' \n']);
fprintf(1, [num2str(min(P_cord(:,3))),' , ',  num2str(max(P_cord(:,3))), ' \n']);
%-----------
n_per_data = 4000; 
n_per_start= 200;
n_per      = 4000;

T_set      = (n_per_start-1)*100+1:(n_per-1)*100+1;   %%%% change time interval %%%%
T_set_dns  = n_per_start:n_per;  

% n_per_data = 1000;
% n_per_start= 1;
% n_per      = 4000;
% T_set      = (n_per_start-1)*100+1:100:(n_per-1)*100+1;

%----------            
for i=1:length(Models)
    %DNS----------------
    if i==1&&cmp_res(i)==1
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_'];
        Cur_name = 'dns_filter04';
        savename = [Dir_save, Cur_name, '_SpAvg', '.mat'];
        str1 = 'b-'; str2=[str2;'DNS        ']; tt=tt+1;
        
        if Load_FFT(i)==1
            load(savename)
        else
            load Matrices/r6/SNP_COEFF.dat
            d  = SNP_COEFF(:, 2:7); a  = d';
            skip = 1;
            %----------------------------------------------
            a = a(:,T_set_dns);
            [N, M]=size(a); time = [1:skip:M]*delta_t/skip;
            fprintf(1,['DNS-The size of a is ', num2str(N), 'by',...
                num2str(M), '\n']);
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'DNS');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
    %EF-POD-----------------
    elseif i==2&&cmp_res(i)==1
%         de=0.001367;
%         Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
%         Cur_name = ['EFROM_new'] ;
%         savename = [Dir_save, Cur_name, num2str(de), '_a_spAvg', '.mat'];
%         str1 = 'r--';  str2=[str2;'EF-ROM-DF  ']; tt=tt+1;
%         if Load_FFT(i)==1
%             load(savename)
%         else
%             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
%                 num2str(n_per),'EFROM_new'];
%             Loadname = [Dir_load, num2str(de), '_a.mat'];
%             load(Loadname, 'a');
%             a = a(:, T_set);
%             [N, M]=size(a); skip = 100;
%             fprintf(1,['EF-The size of a is ', num2str(N), 'by',...
%                 num2str(M),' skip', num2str(skip), '\n']);
%             delta_t = 0.075/skip;
%             %delta_t = 0.01875/skip;
%             time = [1:skip:M]*delta_t/skip;
%             %----------------- FFT
%             [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
%                 Energy_probe, Omega_temp]=FFT_EgySpm(...
%                 a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
%             save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
%                 'Omega_temp', 'Probe_index', 'x_cord');
%         end
        de=0.3;
        nu=0.0285;
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
        Cur_name = ['ADROM_new_FEDF'] ;
        savename = [Dir_save, Cur_name, num2str(de),'_',num2str(nu),'_a_spAvg', '.mat'];
        str1 = 'r--';  str2=[str2;'AD-ROM     ']; tt=tt+1;
        if Load_FFT(i)==1
            load(savename)
        else
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per),'ADROM_new_FEDF'];
            Loadname = [Dir_load, num2str(de),'_',num2str(nu),'_a.mat'];
            load(Loadname, 'a');
            a = a(:, T_set);
            [N, M]=size(a); skip = 100;
            fprintf(1,['ADL-The size of a is ', num2str(N), 'by',...
                num2str(M),' skip', num2str(skip), '\n']);
            delta_t = 0.075/skip;
            time = [1:skip:M]*delta_t/skip;
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
        
     elseif i==3&&cmp_res(i)==1
%         r1=4;
%         Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
%         Cur_name = ['EFProj_'] ;
%         savename = [Dir_save, Cur_name, num2str(r1), '_a_spAvg', '.mat'];
%         str1 = 'r--';  str2=[str2;'EF-ROM-Proj']; tt=tt+1;
%         if Load_FFT(i)==1
%             load(savename)
%         else
%             Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
%                 num2str(n_per),'EFProj_'];
%             Loadname = [Dir_load, num2str(r1), '_a.mat'];
%             load(Loadname, 'a');
%             a = a(:, T_set);
%             [N, M]=size(a); skip = 100;
%             fprintf(1,['EFProj-The size of a is ', num2str(N), 'by',...
%                 num2str(M),' skip', num2str(skip), '\n']);
%             delta_t = 0.075/skip;
%             time = [1:skip:M]*delta_t/skip;
%             %----------------- FFT
%             [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
%                 Energy_probe, Omega_temp]=FFT_EgySpm(...
%                 a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
%             save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
%                 'Omega_temp', 'Probe_index', 'x_cord');
%         end

        de=0.4;
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
        Cur_name = ['ED'] ;
        savename = [Dir_save, Cur_name, num2str(de),'_a_spAvg', '.mat'];
        str1 = 'r--';  str2=[str2;'ED-ROM     ']; tt=tt+1;
        if Load_FFT(i)==1
            load(savename)
        else
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per),'ED'];
            Loadname = [Dir_load, num2str(de),'_a.mat'];
            load(Loadname, 'a');
            a = a(:, T_set);
            [N, M]=size(a); skip = 100;
            fprintf(1,['ADL-The size of a is ', num2str(N), 'by',...
                num2str(M),' skip', num2str(skip), '\n']);
            delta_t = 0.075/skip;
            time = [1:skip:M]*delta_t/skip;
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
        elseif i==4&&cmp_res(i)==1
        de=0.247;
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_1000'];
        Cur_name = ['LDF_'] ;
        savename = [Dir_save, Cur_name, num2str(de), '_a_spAvg', '.mat'];
        str1 = 'r--';  str2=[str2;'L-ROM-DF   ']; tt=tt+1;
        if Load_FFT(i)==1
            load(savename)
        else
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per),'LDF_'];
            Loadname = [Dir_load, num2str(de), '_a.mat'];
            load(Loadname, 'a');
            a = a(:, T_set);
            [N, M]=size(a); skip = 100;
            fprintf(1,['LDF-The size of a is ', num2str(N), 'by',...
                num2str(M),' skip', num2str(skip), '\n']);
            delta_t = 0.075/skip;
            time = [1:skip:M]*delta_t/skip;
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
        elseif i==5&&cmp_res(i)==1
        r1=1;
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
        Cur_name = ['LProj_'] ;
        savename = [Dir_save, Cur_name, num2str(r1), '_a_spAvg', '.mat'];
        str1 = 'r--';  str2=[str2;'L-ROM-Proj ']; tt=tt+1;
        if Load_FFT(i)==1
            load(savename)
        else
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per),'LProj_'];
            Loadname = [Dir_load, num2str(r1), '_a.mat'];
            load(Loadname, 'a');
            a = a(:, T_set);
            [N, M]=size(a); skip = 100;
            fprintf(1,['LProj-The size of a is ', num2str(N), 'by',...
                num2str(M),' skip', num2str(skip), '\n']);
            delta_t = 0.075/skip;
            time = [1:skip:M]*delta_t/skip;
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
        
        elseif i==6&&cmp_res(i)==1
        Dir_save = ['Matrices/energy_spectrum/','energy_spectrum_4000'];
        Cur_name = ['_GROM'] ;
        savename = [Dir_save, Cur_name, '_a_spAvg_ode45', '.mat'];
        str1 = 'r--';  str2=[str2;'G-ROM      ']; tt=tt+1;
        if Load_FFT(i)==1
            load(savename)
        else
            Dir_load = ['Matrices/r',num2str(q_dim),'/connective_matrix145_193_17_',...
                num2str(n_per),'_GROM'];
            Loadname = [Dir_load, '_a.mat'];
            load(Loadname, 'a');
            a = a(:, T_set);
            [N, M]=size(a); skip = 100;
            fprintf(1,['GROM-The size of a is ', num2str(N), 'by',...
                num2str(M),' skip', num2str(skip), '\n']);
            delta_t = 0.075/skip;
            time = [1:skip:M]*delta_t/skip;
            %----------------- FFT
            [ff, Omega, Omega_u, Omega_v, Omega_w, NFFT, time,...
                Energy_probe, Omega_temp]=FFT_EgySpm(...
                a, delta_t, skip, Probe_index, do_avg, do_fig,'POD');
            save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
                'Omega_temp', 'Probe_index', 'x_cord');
        end
     end
    
    %----------------
    % plot energy spectrum
    figure(1)
    set(gcf,'Units','pixels','Position',[1 1000 900 300],'Units','inches')
    set(gcf, 'PaperPosition',[1 2 6 2])
    legend_idx(1,tt)=loglog(ff, 2*abs(Omega(1:NFFT/2+1)), str1,'linewidth',1.1); hold on
    if i==length(Models)
        legend_handle= legend(legend_idx(1,:), str2);
        set(legend_handle,'FontSize',10)
        figure_Egy = ['Per_',num2str(n_per_data),'_',cmp_model,str_avg];
        set(gca,'FontSize',11)
        set(gca,'FontWeight','bold')
        print('-f', '-depsc',figure_Egy)
    end
    xlabel('frequency')
    ylabel('kinetic energy')
end