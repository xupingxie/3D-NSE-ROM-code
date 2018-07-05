clear,clc
% plot energy spectrum for DNS
% if do_avg ==0, compute energy spectrum only on one point
% id do_avg ==1, compute average energy spectrum on several points
% Do FFT on kenectic energy - 0.5*(u^2+v^2+w^2)
load Matrices/r6/connective_matrix145_193_17.mat
%load Matrices/r6/SNP_COEFF.dat
load Matrices/r6/filter_DNS_de04_FE.mat
a=d_filter;
%d  = SNP_COEFF(:, 2:7); a  = d'; 
t_dns=0.075*(1:2000);
Dir_save = 'Matrices/energy_spectrum/energy_spectrum_';
Cur_name = 'dns_filter04'; % change this

delta_t = 0.075; 
skip =1;
do_fig = 0;
do_avg   = 1; % choose do_avg to be 1 and make spatial 
              % avg to draw egy spectrum.
legend_fig = 0;
x_cord = [0.9992, 0.3575, 1.0625]; 

if do_avg==0
    Probe_index = find(abs(x(:,1)-x_cord(1))<0.001&abs(x(:,2)-x_cord(2))<0.001...
        &abs(x(:,3)-x_cord(3))<0.001); % just one point
    savename = [Dir_save, Cur_name, '_', num2str(Probe_index), '.mat'];
else
    Probe_index = find(abs(x(:,1)-x_cord(1))<0.1...
        &abs(x(:,2)-x_cord(2))<0.1...
        &abs(x(:,3)-x_cord(3))<0.2); % choose a box then take avg.
    savename = [Dir_save, Cur_name, '_SpAvg.mat'];
end

P_cord = x(Probe_index,:);
fprintf(1, [num2str(min(P_cord(:,1))),' , ', num2str(max(P_cord(:,1))), ' \n']);
fprintf(1, [num2str(min(P_cord(:,2))),' , ',  num2str(max(P_cord(:,2))), ' \n']);
fprintf(1, [num2str(min(P_cord(:,3))),' , ',  num2str(max(P_cord(:,3))), ' \n']);

%[ff, Omega, NFFT, time, Energy_probe, Omega_temp]=FFT_EgySpm(...
%    a, delta_t, skip, Probe_index, do_avg, do_fig);

[ff, Omega, Omega_u, Omega_v, Omega_w, ...
    NFFT, time, Energy_probe, Omega_temp]=FFT_EgySpm(...
    a, delta_t, skip, Probe_index, do_avg, do_fig, legend_fig);

save(savename, 'ff', 'Omega', 'NFFT', 'time','Energy_probe', 'P_cord', ...
    'Omega_temp', 'Probe_index', 'x_cord');

figure(1)
for i_p = 1: length(Probe_index)
    loglog(ff, 2*abs(Omega_temp(1:NFFT/2+1, i_p))), 
    hold on
    xlabel('frequency')
    ylabel('magnitude')
end
title(num2str(length(Probe_index)))
figure_Egy = ['Energy_figs/DNS', num2str(length(Probe_index))];
%print('-f', '-depsc',figure_Egy)
hold off

figure(2)
loglog(ff, 2*abs(Omega(1:NFFT/2+1)))
xlabel('frequency')
ylabel('magnitude')
title('After average')
figure_Egy = ['Energy_figs/DNS_Avg', num2str(length(Probe_index))];
%print('-f', '-depsc',figure_Egy)