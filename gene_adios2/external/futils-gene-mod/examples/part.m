%
%   Total of steps from attribute "nsteps"
%
nsteps=hdf5read('part.h5','/part/nsteps');
%
%   Particle coordinates at step = nsteps
%
name=sprintf('/part/%3.3d', nsteps);
coor = hdf5read('part.h5',name);

r=coor(:,1);
z=coor(:,2);
phi=coor(:,3);

x=r.*cos(phi);
y=r.*sin(phi);

scatter3(x,y,z,'.'); axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
%
% Get time attribute of dataset
%
attr=sprintf('%s/time', name);
time=hdf5read('part.h5',attr);
title(sprintf('Particles coordinates at time = %.1f',time))
