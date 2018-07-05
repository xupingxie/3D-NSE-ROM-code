function [S11,S12,S13,S22,S23,S33,T] = generate_matrices3d(POD_u, POD_v, POD_w, x, e_conn)

[N         ,r      ] = size(POD_u );
[n_elements,nel_dof] = size(e_conn);
% div_p=0;
max_elem_per_partition = 1000;

n_part = floor( n_elements / (max_elem_per_partition+1) ) + 1;
elem_segment = floor( linspace(0,n_elements,n_part+1) );
max_part_size = max( diff( elem_segment ) );

n_gauss = 5;
[rr,ss,tt,wt] = threed_gauss(n_gauss);
one        = ones(length(rr),1);

for n_pt = 1:n_part
    II   = sparse( max_part_size*nel_dof^2,1 );
    JJ   = II;
    XS11 = sparse( max_part_size*nel_dof^2,1 );
    XS12 = XS11;
    XS13 = XS11;
    %     XS21 = zeros( max_part_size*nel_dof^2,1 );
    XS22 = XS11;
    XS23 = XS11;
    %     XS31 = zeros( max_part_size*nel_dof^2,1 );
    %     XS32 = XS31;
    XS33 = XS11;
    
    for i=1:r
        XT(i).ux = sparse( max_part_size*nel_dof^2,1 );
        XT(i).uy = sparse( max_part_size*nel_dof^2,1 );
        XT(i).uz = sparse( max_part_size*nel_dof^2,1 );
        XT(i).vx = sparse( max_part_size*nel_dof^2,1 );
        XT(i).vy = sparse( max_part_size*nel_dof^2,1 );
        XT(i).vz = sparse( max_part_size*nel_dof^2,1 );
        XT(i).wx = sparse( max_part_size*nel_dof^2,1 );
        XT(i).wy = sparse( max_part_size*nel_dof^2,1 );
        XT(i).wz = sparse( max_part_size*nel_dof^2,1 );
    end
    
    entry_counter = 0;
    
    for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
        nodes_local             = e_conn(n_el,:);
        x_local                 = x(nodes_local,:);
        [x_g,wt_g,phi,phi_x,phi_y,phi_z] = threed_shape(x_local,rr,ss,tt,wt);
        
        S11_loc = threed_bilinear( one, phi_x, phi_x, wt_g) ;
        S12_loc = threed_bilinear( one, phi_x, phi_y, wt_g) ;
        S13_loc = threed_bilinear( one, phi_x, phi_z, wt_g) ;
%               S21_loc = threed_bilinear( one, phi_y, phi_x, wt_g) ;
        S22_loc = threed_bilinear( one, phi_y, phi_y, wt_g) ;
        S23_loc = threed_bilinear( one, phi_y, phi_z, wt_g) ;
        %       S31_loc = threed_bilinear( one, phi_z, phi_x, wt_g) ;
        %       S32_loc = threed_bilinear( one, phi_z, phi_y, wt_g) ;
        S33_loc = threed_bilinear( one, phi_z, phi_z, wt_g) ;
        
%         podu_xk= phi_x*POD_u(nodes_local,:); podv_yk=
%         phi_y*POD_v(nodes_local,:);
%         div_p  = div_p + twod_bilinear(one, podu_xk+podv_yk, one, wt_g);
        
        %  Assemble contributions into the global system matrix
        %-------------------------------------------------------------------
        for i=1:nel_dof
            for j=1:nel_dof
                entry_counter = entry_counter + 1;
                II  (entry_counter) = nodes_local(i);
                JJ  (entry_counter) = nodes_local(j);
                XS11(entry_counter) = S11_loc(i,j);
                XS12(entry_counter) = S12_loc(i,j);
                XS13(entry_counter) = S13_loc(i,j);
                XS22(entry_counter) = S22_loc(i,j);
                XS23(entry_counter) = S23_loc(i,j);
                XS33(entry_counter) = S33_loc(i,j);
                
                for k=1:r   % recomputing the local arrays extra times, should break out
                    % this part of the calculation...
                    podu_k = phi*POD_u(nodes_local,k);
                    podv_k = phi*POD_v(nodes_local,k);
                    podw_k = phi*POD_w(nodes_local,k);
                    
                    Tux_local = threed_bilinear( podu_k, phi_x, phi, wt_g );
                    Tuy_local = threed_bilinear( podu_k, phi_y, phi, wt_g );
                    Tuz_local = threed_bilinear( podu_k, phi_z, phi, wt_g );
                    Tvx_local = threed_bilinear( podv_k, phi_x, phi, wt_g );
                    Tvy_local = threed_bilinear( podv_k, phi_y, phi, wt_g );
                    Tvz_local = threed_bilinear( podv_k, phi_z, phi, wt_g );
                    Twx_local = threed_bilinear( podw_k, phi_x, phi, wt_g );
                    Twy_local = threed_bilinear( podw_k, phi_y, phi, wt_g );
                    Twz_local = threed_bilinear( podw_k, phi_z, phi, wt_g );
                    
                    
                    XT(k).ux(entry_counter) = Tux_local(i,j);
                    XT(k).uy(entry_counter) = Tuy_local(i,j);
                    XT(k).uz(entry_counter) = Tuz_local(i,j);
                    XT(k).vx(entry_counter) = Tvx_local(i,j);
                    XT(k).vy(entry_counter) = Tvy_local(i,j);
                    XT(k).vz(entry_counter) = Tvz_local(i,j);
                    XT(k).wx(entry_counter) = Twx_local(i,j);
                    XT(k).wy(entry_counter) = Twy_local(i,j);
                    XT(k).wz(entry_counter) = Twz_local(i,j);
                end
            end
        end
        if mod(n_el, 50) ==0 && mod(n_el, 500)~=0
            fprintf(1, [num2str(n_el), '\t'])
        elseif mod(n_el, 500) == 0
            fprintf(1, [num2str(n_el), '\n'])
        end
        
    end
    
    if ( n_pt==1 )
%         S11 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS11(1:entry_counter),...
%             N, N );
%         S12 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS12(1:entry_counter),...
%             N, N );
%         S13 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS13(1:entry_counter),...
%             N, N );
%         S22 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS22(1:entry_counter),...
%             N, N );
%         S23 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS23(1:entry_counter),...
%             N, N );
%         S33 = sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS33(1:entry_counter),...
%             N, N );
        
        for k=1:r
            T(k).ux = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).ux(1:entry_counter),...
                N, N );
            T(k).uy = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).uy(1:entry_counter),...
                N, N );
            T(k).uz = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).uz(1:entry_counter),...
                N, N );
            T(k).vx = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vx(1:entry_counter),...
                N, N );
            T(k).vy = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vy(1:entry_counter),...
                N, N );
            T(k).vz = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vz(1:entry_counter),...
                N, N );
            T(k).wx = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wx(1:entry_counter),...
                N, N );
            T(k).wy = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wy(1:entry_counter),...
                N, N );
            T(k).wz = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wz(1:entry_counter),...
                N, N );
        end
        
    else
        
%         S11 = S11 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS11(1:entry_counter),...
%             N, N );
%         S12 = S12 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS12(1:entry_counter),...
%             N, N );
%         S13 = S13 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS13(1:entry_counter),...
%             N, N );
%         S22 = S22 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS22(1:entry_counter),...
%             N, N );
%         S23 = S23 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS23(1:entry_counter),...
%             N, N );
%         S33 = S33 + sparse( II(1:entry_counter), JJ(1:entry_counter),...
%             XS33(1:entry_counter),...
%             N, N );
        
        for k=1:r
            T(k).ux = T(k).ux + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).ux(1:entry_counter),...
                N, N );
            T(k).uy = T(k).uy + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).uy(1:entry_counter),...
                N, N );
            T(k).uz = T(k).uz + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).uz(1:entry_counter),...
                N, N );
            T(k).vx = T(k).vx + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vx(1:entry_counter),...
                N, N );
            T(k).vy = T(k).vy + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vy(1:entry_counter),...
                N, N );
            T(k).vz = T(k).vz + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).vz(1:entry_counter),...
                N, N );
            T(k).wx = T(k).wx + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wx(1:entry_counter),...
                N, N );
            T(k).wy = T(k).wy + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wy(1:entry_counter),...
                N, N );
            T(k).wz = T(k).wz + sparse( II(1:entry_counter), JJ(1:entry_counter),...
                XT(k).wz(1:entry_counter),...
                N, N );
        end
    end
    
end

% end  generate_matrices3d

function entry=iszero(entry)
if abs(entry)<eps
    entry =0;
end