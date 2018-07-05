function [M] = compute_3Dmass_matrix(x,e_conn)
%================================================================================
% compute mass matrices for the 3D case
%================================================================================
[N,            n_dim] = size(x);
[n_elements, nel_dof] = size(e_conn);

max_elem_per_partition = 10000;

n_part = floor( n_elements / (max_elem_per_partition+1) ) + 1;
elem_segment = floor( linspace(0,n_elements,n_part+1) );
max_part_size = max( diff( elem_segment ) );

n_gauss = 5;
[rr,ss,tt,wt] = threed_gauss(n_gauss);
one        = ones(length(rr),1);

for n_pt=1:n_part
    II        = zeros(nel_dof^2*max_part_size,1);
    JJ        = zeros(nel_dof^2*max_part_size,1);
    XX        = zeros(nel_dof^2*max_part_size,1);
    entry_counter = 0;
    
    for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
        % compute value of each test function and their spatial derivaties
        % at the integration points
        nodes_local                     = e_conn(n_el,:);
        x_local                         = x(nodes_local,:);
        [x_g,wt_g,phi] = threed_shape(x_local,rr,ss,tt,wt);
        
        % compute the value of each variable at the Gauss points
        
        M_loc   = threed_bilinear( one, phi  , phi  , wt_g) ;
        
        %-------------------------------------------------------------------------
        %  Assemble contributions into the global system matrix
        %-------------------------------------------------------------------------
        for i=1:nel_dof
            for j=1:nel_dof
                entry_counter = entry_counter + 1;
                II(entry_counter) = nodes_local(i);
                JJ(entry_counter) = nodes_local(j);
                XX(entry_counter) = M_loc(i,j);
            end
        end
    end
    
    if ( n_pt==1 )
        M = sparse( II(1:entry_counter), JJ(1:entry_counter),...
            XX(1:entry_counter),...
            N, N );
    else
        M = M + ...
            sparse( II(1:entry_counter), JJ(1:entry_counter),...
            XX(1:entry_counter),...
            N, N );
    end
end