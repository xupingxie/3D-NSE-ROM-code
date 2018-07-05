%===========================================================================
% 3D flow past a circular cylinder;
% machine learning PODROM using the ODE 45 solver.
%
%    @Xuping Xie  
% Oak Ridge National Lab
%    10/06/2017  
%===========================================================================
global q_dim n_x n_y n_z

generate_s = 0;     % generate stiffness matrix, etc
generate_t = 0;     % generate matrices for nonlinear term
q_dim      = 6;     % number of POD basis
%------ set several parameters
n_x        = 145;   % number of nodes in x
n_y        = 193;   % number of nodes in y
n_z        = 17 ;   % number of nodes in z
delta_t    = 0.075; % time step size
n_per      = 4000;   % number of steps
fprintf(1, ['compute the first ', num2str(n_per*0.075), ' seconds\n']);
t_range    = 62.5:0.075:62.5+(n_per-1)*0.075;
mu         = 0.001; % 1/Re
load_rom_mat = 1;  % 0= no ROM matrix, need compute

fprintf('Data preprocessing...\n')
% R_c: coarse rate in both x and y direction; R_z: coarse rate in z direction
R_c = 1; R_z = 1;
n_nodesx = (n_x-1)/R_c+1;
n_nodesy = (n_y-1)/R_c+1;
n_nodesz = (n_z-1)/R_z+1;


filename = ['Matrices/r',num2str(q_dim),'/connective_matrix',...
    num2str(n_nodesx),'_',num2str(n_nodesy),'_', num2str(n_nodesz)];
filename1 = [filename, '.mat'];   % file contains POD basis, mesh information
filename2 = [filename, '_S.mat']; % file contains mass ans stiffness matrices

filenamet = [filename, '_T.mat' ];
filenametu= [filename, '_Tu.mat'];
filenametv= [filename, '_Tv.mat'];
filenametw= [filename, '_Tw.mat'];


%------ load initial condition
if q_dim == 4
    load Matrices/r6/SNP_COEFF.dat
    a0 = SNP_COEFF(1, 2:5); a0 = a0';
    d  = SNP_COEFF(:, 2:5); d  = d' ;
elseif q_dim == 6
    load Matrices/r6/SNP_COEFF.dat
    a0 = SNP_COEFF(1, 2:7); a0 = a0';
    d  = SNP_COEFF(:, 2:7); d  = d' ;
end
%---------compute derivate of coeff
filename_rommat = [filename,'GROM_ABC.mat'];
load(filename_rommat)
  
adot=compute_derivative_coeff(d,delta_t);
train=Training_data_compute(adot,d,q_dim,A,B,C);





%generate tecplot files
% load(filename1)
% POD_u = POD_u'; center_u = center_u';
% POD_v = POD_v'; center_v = center_v';
% POD_w = POD_w'; center_w = center_w';
% for fig = 1:100;
%     u1 = center_u + POD_u*a(:,fig);
%     v1 = center_v + POD_v*a(:,fig);
%     w1 = center_w + POD_w*a(:,fig);
%     plt_name = ['tecplot_figures/tec', num2str(fig,'%05d'), '.plt'];
%     threed_to_tecplot(x,e_conn,[u1,v1,w1],plt_name,num2str(fig,'%05d'))
% end

%end
