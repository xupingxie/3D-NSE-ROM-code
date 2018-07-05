function [ff, Omega, Omega_u, Omega_v, Omega_w, ...
    NFFT, time, Energy_probe, Omega_temp]=FFT_EgySpm(...
    a, delta_t, skip, Probe_index, do_avg, do_fig, legend_fig)
% plot energy spectrum for DNS
% Do FFT on kenectic energy - 0.5*(u^2+v^2+w^2)
[N, M]=size(a); 

load(['Matrices/r',num2str(N),'/connective_matrix145_193_17.mat'])

time = [1:skip:M]*delta_t;
delta_t=skip*delta_t;

POD_u = POD_u'; center_u = center_u';
POD_v = POD_v'; center_v = center_v';
POD_w = POD_w'; center_w = center_w';
t_j=1;
for t_i = 1:skip:M;
    u = center_u + POD_u*a(:,t_i);
    v = center_v + POD_v*a(:,t_i);
    w = center_w + POD_w*a(:,t_i);
    U_probe(t_j).mat = [u(Probe_index,1)'; v(Probe_index,1)'; w(Probe_index,1)'];
    if do_avg==0
        Energy_probe(t_j) = sum(0.5*U_probe(t_j).mat.*U_probe(t_j).mat);
        u_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(1,:).*U_probe(t_j).mat(1,:);
        v_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(2,:).*U_probe(t_j).mat(2,:);
        w_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(3,:).*U_probe(t_j).mat(3,:);
    else
        Energy_probe(t_j,:) = sum(0.5*U_probe(t_j).mat.*U_probe(t_j).mat);
        u_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(1,:).*U_probe(t_j).mat(1,:);
        v_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(2,:).*U_probe(t_j).mat(2,:);
        w_Egy_probe(t_j, :) = 0.5*U_probe(t_j).mat(3,:).*U_probe(t_j).mat(3,:);
    end
    
    t_j=t_j+1;
    if mod(t_i-1,1000)==0
        fprintf(1, [num2str(t_i), '\n'])
    end
end

NFFT = 2^nextpow2(length(Energy_probe));
if do_avg==0
    Omega=fft(Energy_probe, NFFT)/length(Energy_probe);
    Omega_u=fft(u_Egy_probe, NFFT)/length(u_Egy_probe);
    Omega_v=fft(v_Egy_probe, NFFT)/length(v_Egy_probe);
    Omega_w=fft(w_Egy_probe, NFFT)/length(w_Egy_probe);
    Omega_temp=Omega;
else
    for i_p=1:length(Probe_index)
        Omega_temp(:, i_p)= fft((Energy_probe(:, i_p)), NFFT)/length(Energy_probe); % take square root?
        Omega_u_tmp(:, i_p)= fft((u_Egy_probe(:, i_p)), NFFT)/length(u_Egy_probe);
        Omega_v_tmp(:, i_p)= fft((v_Egy_probe(:, i_p)), NFFT)/length(v_Egy_probe);
        Omega_w_tmp(:, i_p)= fft((w_Egy_probe(:, i_p)), NFFT)/length(w_Egy_probe);
    end
    Omega=sum(Omega_temp')'/length(Probe_index);
    Omega_u=sum(Omega_u_tmp')'/length(Probe_index);
    Omega_v=sum(Omega_v_tmp')'/length(Probe_index);
    Omega_w=sum(Omega_w_tmp')'/length(Probe_index);
end
ff = (1/delta_t)/2*linspace(0,1,NFFT/2+1);

if do_fig==1
    figure
    % loglog(ff(1:length(plot_half)), spectrum(1:length(plot_half)) )
    loglog(ff, 2*abs(Omega(1:NFFT/2+1)), 'b--')
    xlabel('frequency')
    ylabel('magnitude')
    legend(legend_fig)
    title('After average')
end

