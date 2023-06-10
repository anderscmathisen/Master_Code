clear all
home

res = '-r150';
plotx2east
plotzIntoPlane

dir_sample = 'Datapath';
dir_mtex = 'Datapath/plots';

lattice_parameters = [6.1 6.1 11];
angles = [90 90 120] * degree;

cs = crystalSymmetry('6/mmm', lattice_parameters, angles, 'mineral',...
    'ErMnO3', 'X||a*', 'Z||c*'); % Note the X||a* alignment



directions = {xvector, yvector, zvector};

ebsd = EBSD.load(fullfile(dir_sample,...
        "ErMnO3_full_lamella_rotated.ang"), cs,...
         'setting 2');

rot = rotation.byAxisAngle(zvector, 180*degree)
ebsd = rotate(ebsd, rot, 'keepXY')


ebsd.scanUnit = 'nm';
ebsd(ebsd.ci < 0.001).phase = -1;

ipfkey = ipfTSLKey(ebsd.CS);



%reconstruct grains
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle',...
        10*degree);
grains_ermno3 = grains('ErMnO3');

big_grains = grains_ermno3(grains_ermno3.grainSize > 130); %exclude smallest grains

ipfkey.inversePoleFigureDirection = directions{3};
omcolor = ipfkey.orientation2color(ebsd.orientations);

cS = crystalShape.hex(ebsd.CS);


%plot figure labeling grains, and print to console mean orientation of grains
% results are copy pasted to .txt file and read by Gonio_tilt.ipynb
% Results are in Full_lamella_oris_rotated180.txt, with renamed labels according to thesis 
% and deleted uninteresting grains. 

figure
plot(ebsd, omcolor)

hold on
for i=1 : length(big_grains) 
    try
        text(big_grains(i),(i), 'fontsize' ,16)
        ori = big_grains(i);
        fprintf("%i: [%d,%d, %d, %d]\n", i, ori.meanOrientation.a, ori.meanOrientation.b, ori.meanOrientation.c, ori.meanOrientation.d)
    catch exception
    end
end
legend('off')