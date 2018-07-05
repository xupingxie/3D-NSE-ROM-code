function [ element, iso1, iso2 ] = mesh_search( barycenter, xy, e_conn, point )
%%
%  MESH_SEARCH  A MATLAB function that locates which element contains a
%  given point.  The function only provides one such element and makes an
%  explicit assumption that the triangulation is Delauney.
%%
found = 0;
radii = [0.5];
index_search = 1:size(e_conn,1);
index_elem_0 = [];
for k=1:length(radii)
    index_elem = index_search(find(abs(barycenter(index_search,1)-point(1))<radii(k) &...
        abs(barycenter(index_search,2)-point(2))<radii(k)));
    if isempty(index_elem)==1
        element = 0; iso1= 0; iso2= 0;
        break
    elseif length(index_elem)==1
        break
    else
        index_search = index_elem;
        index_elem_0 = index_elem;
    end
end
if isempty(index_elem_0)~=1 && isempty(index_elem)==1
    index_elem = index_elem_0;
end
% [n_elem,npe] = size(e_conn);
if isempty(index_elem)==1
    element = 0; iso1= 0; iso2= 0;
else
    for i = 1:length(index_elem)
        element = index_elem(i);
        x_local = xy( e_conn(element,1:3), : );
        
        %    plot( [ x_local(1,1) x_local(2,1) ],[ x_local(1,2) x_local(2,2)],...
        %          [ x_local(2,1) x_local(3,1) ],[ x_local(2,2) x_local(3,2)],...
        %          [ x_local(3,1) x_local(1,1) ],[ x_local(3,2) x_local(1,2)] )
        
        delta = ( x_local(3,1)-x_local(2,1) )*( x_local(1,2)-x_local(2,2) ) ...
            - ( x_local(3,2)-x_local(2,2) )*( x_local(1,1)-x_local(2,1) );
        
        iso1  = ( ( x_local(2,1)-point(1) )*( x_local(3,2)-point(2) ) ...
            - ( x_local(3,1)-point(1) )*( x_local(2,2)-point(2) ) ) / delta;
        
        iso2  = ( ( x_local(3,1)-point(1) )*( x_local(1,2)-point(2) ) ...
            - ( x_local(1,1)-point(1) )*( x_local(3,2)-point(2) ) ) / delta;
        
        iso3  = ( ( x_local(1,1)-point(1) )*( x_local(2,2)-point(2) ) ...
            - ( x_local(2,1)-point(1) )*( x_local(1,2)-point(2) ) ) / delta;
        
        if ( iso1>0 & iso2>0 & iso3>0 )
            found = 1;
            return
        end
    end
    
end
if found == 0
    element = 0; iso1= 0; iso2= 0;
end

end