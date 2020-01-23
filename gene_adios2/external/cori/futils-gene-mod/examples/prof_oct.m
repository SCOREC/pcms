#!/usr/bin/octave
load "prof.h5" "profile_1d"

nsteps=size(profile_1d.time,2);
hold on
for i=1:nsteps
   plot(profile_1d.xgrid,profile_1d.phi(:,i));
end
xlabel('X');
ylabel('\phi');
title('Time Evolution');
