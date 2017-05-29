cmap = zeros(41,3);
cmap(:,1) = [0:.05:1 ones(1,20)]';       % R   fades up to white, then solid
cmap(:,2) = [0:.05:1 .95:-.05:0]';      % G   fades jup to white then fadces down
cmap(:,3) = [ones(1,20) 1:-.05:0]'      % B   solid then fades down