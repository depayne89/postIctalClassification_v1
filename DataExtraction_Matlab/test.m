
sz_lengths = [ 1 2 3 4 26 27 29 31];    %make mock sz lengths
sz_num = 1:size(sz_lengths,2);          % make x-axis values (sz number)
figure

subplot(1,2,1)                          % plot to left of figure
scatter(sz_lengths, sz_num, 8, 'filled')             % plot sz_lengths
axis([0 max(sz_lengths)*1.1 .5 size(sz_lengths,2)+.5]); % adjust axis to fit best
set(gca, 'ycolor', 'none', 'ydir', 'reverse', 'position', [0.04 0.1 0.15 0.8]); % ycolor removes y axes, yir and xdir reverses axes, position sets position in [left, bottom, width, height]
times = 10;
rand_sz = repmat(1:times,[size(sz_lengths,2) 1]) + rand(size(sz_lengths,2), times);     % make mock data
xlabel('Seizure length', 'Fontname', 'Calibri', 'Fontsize', 8);
set(gca, 'FontName', 'Calibri', 'Fontsize', 6)

subplot(1,2,2)      % allocate to right of image
imagesc(rand_sz)    % plot survival curve
                               

set(gca, 'ycolor', 'none', 'box', 'off', 'position', [0.21 0.1 0.68 0.8]); % y color to none removes axis, box off removes top x axis, position sets position in [left, bottom, width, height]

xlabel('Time since seizure end (s)', 'Fontname', 'Calibri', 'Fontsize', 8);
set(gca, 'FontName', 'Calibri', 'Fontsize', 6)
cmap = zeros(41,3);
% Start: [ 0/256 34/256 88/256   Middle [1 1 1] End  [164/256 3/256 53/256]
% 
cmap(:,1) = [0:      1/20:         1   1-92/(20*256):  -92/(20*256):  164/256]';          % Red   fades up to white, then solid
cmap(:,2) = [34/256: 222/(20*256): 1   1-253/(20*256): -253/(20*256): 3/256]';          % Green  fades jup to white then fadces down
cmap(:,3) = [88/256: 168/(20*256): 1   1-203/(20*256): -203/(20*256): 53/256]';      	% Blue   solid then fades down
colormap(cmap);                             % sets colormap to above
colorbar('Position', [.9 0.096 .05 .76]);    % shows colormap

h = suptitle('Patient 11')
set(h,'FontSize',8,'Fontname','calibri')

fig=gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];
% print -dtiff test600.tif -r600
print(['Figures/ScatSurv/surv_' num2str(11)], '-dtiff', '-r600')   % -rx changes dpi/ppi to x

