%
%     Get data from data sets
%
time=hdf5read('prof.h5','/profile_1d/time');
x=hdf5read('prof.h5','/profile_1d/xgrid');
phi=hdf5read('prof.h5','/profile_1d/phi');
%
%     Attributes
%
attr=hdf5read('prof.h5','/title'); title_ann=attr.Data;
attr=hdf5read('prof.h5','/profile_1d/time/title'); time_ann=attr.Data;
attr=hdf5read('prof.h5','/profile_1d/xgrid/title'); x_ann=attr.Data;
attr=hdf5read('prof.h5','/profile_1d/phi/title'); phi_ann=attr.Data;

figure
nsteps=size(time,2)
hold on
for i = 1:10
    plot(x,phi(i,:));
end
xlabel(x_ann)
ylabel(phi_ann)
title(title_ann)

figure
[c,h]=contourf(double(x), double(time(1:10)), double(phi(1:10,:))); clabel(c,h); colorbar
xlabel(time_ann)
ylabel(x_ann)
title(phi_ann)
