function [Dweight1, Dweight2, Origin1, Origin2, Terminus1, Terminus2, ind_nz,boundary] ...
          = set_interfaces(Slice)

%Due to the many complexities of dealing with a brain embedded in a
%square... Use this function to determine the important pixels and keep
%record of them. 

%Dweights appear to hold an average pixel intensity value, I believe to
%account for the different diffusion parameters in gray and white matter
%(ahd 6-30-2011)

%Origin{1,2} and Terminus{1,2} seem to hold the proper indicies to
%efficiently caluclate matrix values. 1 corresponds to columns (y
%direction), 2 to rows (x direction)

%ind_nz = all pixels to be considered as nodes

%boundary = nodes that are boundary nodes

% Read the dimensions of the brain slice
    n1 = size(Slice, 1);
    n2 = size(Slice, 2);

% Get the indices of cells with nonzero intensity,
% and find how many of these there are.
   
    ind_nz = find(Slice > 0);
    num_nz = length(ind_nz);

    %counter for boundary nodes
    b = 0;


% Set up interfaces along the first dimension.
    Origin1 = 0;
    Terminus1 = 0;
    Dweight1 = 0;
    j=0;
    for i = 1:num_nz
      I = ind_nz(i);
      %ahd-7-1-2011
      %seems that this check is making sure to not include boundaries?
      %First check is to guarantee not on the bottom border of the image. Second check 
      %is to make sure pixel below current one is "in brain." 
      %(Why no check for top border of for above pixel?)Aha! Think of
      %center of pixel as node, and it makes sense.
      %Maybe relevant for Boundary conditions??? 
      if (mod(I,n1) ~= 0) && ismember(I+1,ind_nz)
        j = j+1;
        Origin1(j,1) = I;
        Terminus1(j,1) = I+1;
        Dweight1(j,1) = (Slice(I) + Slice(I+1))/2;
        %Lweight1(j,1) = (CSFSlice(I) + CSFSlice(I+1))/2;
      end
    end

    
% Set up interfaces along the second dimension
    Origin2 = 0;
    Terminus2 = 0;
    Dweight2 = 0;
    j=0;
    for i = 1:num_nz
      I = ind_nz(i);
      %ahd-7-1-2011 
      %I would think this should be similar check as above... First check would just
      %ensure that it isn't in the last column? and second check that the
      %pixel to the right is "in brain." Again... why not check pixel to
      %left or first column?
      if (I+n1 < n1*n2) && ismember(I+n1,ind_nz)
        j=j+1;
        Origin2(j,1) = I;
        Terminus2(j,1) = I+n1;
        Dweight2(j,1) = (Slice(I) + Slice(I+n1))/2;
        %Lweight2(j,1) = (CSFSlice(I) + CSFSlice(I+n1))/2;
      end
    end
    
    %Find Boundary nodes
    for i=1:num_nz
        I = ind_nz(i);
        if ~ismember(I-1,ind_nz) || ~ismember(I+1,ind_nz) || ~ismember(I-n1,ind_nz) || ~ismember(I+n1,ind_nz)
            b = b+1;
            boundary(b) = I;
        end
    end
    