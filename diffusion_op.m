function b = diffusion_op(u, D, dx, ...
                          Weight1, Weight2, ...
                          Origin1, Origin2, ...
                          Terminus1, Terminus2)
                        

Flux1 = D * Weight1 .* (u(Origin1) - u(Terminus1)) / dx^2;
  
Flux2 = D * Weight2 .* (u(Origin2) - u(Terminus2)) / dx^2;

b = zeros(size(u));

% Account for flux leaving along the first dimension
  b(Origin1) = b(Origin1) - Flux1;

  
% Next fluxes entering along the first dimension
  b(Terminus1) = b(Terminus1) + Flux1;
  
% Now flux leaving along the second dimension
  b(Origin2) = b(Origin2) - Flux2;
  
% And flux entering along the second dimension
  b(Terminus2) = b(Terminus2) + Flux2;
