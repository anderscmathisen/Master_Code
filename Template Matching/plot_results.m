clear all
home

res = '-r150';
plotx2east
plotzIntoPlane

dir_sample = '/Users/anders/Library/CloudStorage/OneDrive-NTNU/Master/Data/2100F/27_04_23';
dir_mtex = '/Users/anders/Library/CloudStorage/OneDrive-NTNU/Master/Data/2100F/27_04_23/plots';

lattice_parameters = [6.1 6.1 11];
angles = [90 90 120] * degree;

cs = crystalSymmetry('6/mmm', lattice_parameters, angles, 'mineral',...
    'ErMnO3', 'X||a*', 'Z||c*'); %inverted definition because of bug in Orix



directions = {xvector, yvector, zvector};

ebsd = EBSD.load(fullfile(dir_sample,...
        "ErMnO3_full_lamella_rotated.ang"), cs,...
         'setting 2');

rot = rotation.byAxisAngle(zvector, 180*degree)
ebsd = rotate(ebsd, rot, 'keepXY')


ebsd.scanUnit = 'nm';
ebsd(ebsd.ci < 0.001).phase = -1;

ipfkey = ipfTSLKey(ebsd.CS);

figure
plot(ebsd, ebsd.ci)
mtexColorbar('title', 'NCC')
mtexColorMap black2white
%export_fig(fullfile(dir_mtex, 'maps_ncc.png'), res)


titles = {'x', 'y', 'z'};
for i=1:3
    ipfkey.inversePoleFigureDirection = directions{i};
    omcolor = ipfkey.orientation2color(ebsd.orientations);

    % IPF maps
    figure
    plot(ebsd, omcolor)
    hold on
    %export_fig(fullfile(dir_mtex, ['maps_ipf_' titles{i} '.png']), res)
end

%calculate grains
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle',...
        10*degree);
grains_ermno3 = grains('ErMnO3');

big_grains = grains_ermno3(grains_ermno3.grainSize > 130);

ipfkey.inversePoleFigureDirection = directions{3};
omcolor = ipfkey.orientation2color(ebsd.orientations);

cS = crystalShape.hex(ebsd.CS);

for i=1:3
    
    ipfkey.inversePoleFigureDirection = directions{i};
    omcolor = ipfkey.orientation2color(ebsd.orientations);
    
    figure
    plot(ebsd, omcolor)
    hold on
    %plot(grains_ermno3.boundary)
    hold on
    plot(big_grains, 0.5*cS, 'linewidth', 2, 'colored')
    legend('off')
    %export_fig(fullfile(dir_mtex, ['ipf_with_unitcell_' titles{i} '.png']), res)
    
end


figure
plot(ebsd('indexed'), ebsd('indexed').mis2mean.angle / degree)
caxis([0,20])
mtexColorbar('title', 'Misorientation to mean orientation [deg]')
mtexColorMap LaboTeX
hold on
plot(big_grains.boundary)
