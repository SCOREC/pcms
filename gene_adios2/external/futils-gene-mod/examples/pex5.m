figure;

coor=hdf5read('part_000.h5', 'part');
r=coor(:,1);
z=coor(:,2);
phi=coor(:,3);
x=r.*cos(phi);
y=r.*sin(phi);
subplot(221)
scatter3(x,y,z,'.'); axis equal;xlabel('X');ylabel('Y');zlabel('Z')
title('part_000.h5')

coor=hdf5read('part_001.h5', 'part');
r=coor(:,1);
z=coor(:,2);
phi=coor(:,3);
x=r.*cos(phi);
y=r.*sin(phi);
subplot(222)
scatter3(x,y,z,'.'); axis equal;xlabel('X');ylabel('Y');zlabel('Z')
title('part_001.h5')

coor=hdf5read('part_002.h5', 'part');
r=coor(:,1);
z=coor(:,2);
phi=coor(:,3);
x=r.*cos(phi);
y=r.*sin(phi);
subplot(223)
scatter3(x,y,z,'.'); axis equal;xlabel('X');ylabel('Y');zlabel('Z')
title('part_002.h5')

coor=hdf5read('part_003.h5', 'part');
r=coor(:,1);
z=coor(:,2);
phi=coor(:,3);
x=r.*cos(phi);
y=r.*sin(phi);
subplot(224)
scatter3(x,y,z,'.'); axis equal;xlabel('X');ylabel('Y');zlabel('Z')
title('part_003.h5')
