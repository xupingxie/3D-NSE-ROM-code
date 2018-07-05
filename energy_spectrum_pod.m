clear,clc
% plot energy spectrum for POD-ROM
load Matrices/r6/connective_matrix145_193_17.mat
%Dir_load = 'Matrices/connective_matrix145_193_17_';
Dir_load = 'Matrices/r6/connective_matrix145_193_17_';
Dir_save = 'Matrices/energy_spectrum/energy_spectrum_';
do_avg   = 1; % choose do_avg to be 1 and make spatial avg to draw egy spectrum.
Cur_name = '4000LDF_0.06_a'; % change this
loadname = [Dir_load, Cur_name, '.mat'];
load(loadname);
% load Matrices/connective_matrix145_193_17_2000_a.mat
delta_t=delta_t/100;

if do_avg==0
    Probe_index = find(abs(x(:,1)-1)<0.001&abs(x(:,2)-0.3575)<0.001...
        &abs(x(:,3)-1.0625)<0.001); % just one point
    savename = [Dir_save, Cur_name, '.mat'];
else
    Probe_index = find(abs(x(:,1)-1)<0.001...
        &x(:,2)-1<0.001&x(:,2)+1>0.001...
        &abs(x(:,3)-1.0625)<0.001); % choose several point then take avg.
    savename = [Dir_save, Cur_name, '_SpAvg.mat'];
end



[N, M]=size(a);     skip=100;
T_whole= delta_t*M;
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
    
%     du = center_u + POD_u*d(:,t_i);
%     dv = center_v + POD_v*d(:,t_i);
%     dw = center_w + POD_w*d(:,t_i);
    
    U_probe(t_j,:) = [u(Probe_index,1)', v(Probe_index,1)', w(Probe_index,1)'];
%    DNS_U_probe(t_j,:) = [du(Probe_index,1)', dv(Probe_index,1)', dw(Probe_index,1)'];
    
    if do_avg==0
        Energy_probe(t_j) = 0.5*U_probe(t_j,:)*U_probe(t_j,:)';
    else
        Energy_probe(t_j) = 0.5*U_probe(t_j,:)*U_probe(t_j,:)'/length(Probe_index);
%        DNS_Energy_probe(t_j) = 0.5*DNS_U_probe(t_j,:)*DNS_U_probe(t_j,:)'/length(Probe_index);
    end
    t_j=t_j+1;
    if mod(t_i-1,1000)==0
        fprintf(1, [num2str(t_i), '\n']);
    end
end


NFFT = 2^nextpow2(length(Energy_probe));
Omega=fft(Energy_probe, NFFT)/length(Energy_probe);
ff = (1/delta_t)/2*linspace(0,1,NFFT/2+1);

loglog(ff, 2*abs(Omega(1:NFFT/2+1)), 'r--')
%hold on

%load Matrices/energy_spectrum/energy_spectrum_dns_SpAvg.mat

%dns_NFFT = NFFT;
%dns_Omega= Omega;
%dns_ff = ff;

%loglog(dns_ff,2*abs(dns_Omega(1:dns_NFFT/2+1)), 'b')
%hold off
xlabel('frequency')
ylabel('magnitude')
legend('ROM')
%set(gca,'Position',[0,0,4,2])
save(savename, 'ff', 'Omega', 'NFFT',...
    'time','Energy_probe');

