%
%     Get data from data sets
%
time=hdf5read('sim.h5','/scalars/time');
epot=hdf5read('sim.h5','/scalars/epot');
ekin=hdf5read('sim.h5','/scalars/ekin');
%
%     Attributes
%
attr=hdf5read('sim.h5','/title'); title_ann=attr.Data;
attr=hdf5read('sim.h5','/scalars/time/title'); time_ann=attr.Data;
attr=hdf5read('sim.h5','/scalars/epot/title'); epot_ann=attr.Data;
attr=hdf5read('sim.h5','/scalars/ekin/title'); ekin_ann=attr.Data;

figure
plot(time, ekin, time, epot)
legend(epot_ann, ekin_ann)
xlabel(time_ann)
title(title_ann)
grid on
