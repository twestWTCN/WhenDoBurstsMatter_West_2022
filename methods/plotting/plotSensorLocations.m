function plotSensorLocations(locs,roi,indexSet,powDsync,betaSNR)

% transform to scale
powDsync = -powDsync;
powDsync = (powDsync-min(powDsync)+50)
betaSNR = betaSNR(:,2).^1.4
ofs = 20

% Set 1
subplot(1,3,1)
scatter3(locs(:,1),locs(:,2),locs(:,3),'b')
hold on
scatter3(locs(indexSet{2},1),locs(indexSet{2},2),locs(indexSet{2},3),'b','filled')
scatter3(roi(2,1),roi(2,2),roi(2,3),150,'r','filled','MarkerFaceAlpha',0.8)

plot3([-70 -70+ofs],[60 60],[0 0],'k','LineWidth',2); hold on
plot3([-70 -70],[60 60+ofs],[0 0],'k','LineWidth',2)
plot3([-70 -70],[60 60],[0 0+ofs],'k','LineWidth',2)

view([-125 30]); 
xlim([-50 100]); ylim([-30 50]); zlim([0 80])
axis equal; box off; axis off

% Set 2
subplot(1,3,2)
% scatter3(locs(:,1),locs(:,2),locs(:,3),'k'); 
hold on
scatter3(roi(2,1),roi(2,2),roi(2,3),100,'r','filled')

scatter3(locs(indexSet{2},1),locs(indexSet{2},2),locs(indexSet{2},3),betaSNR,'b')
[~,ia] = intersect(indexSet{2},indexSet{3});
scatter3(locs(indexSet{3},1),locs(indexSet{3},2),locs(indexSet{3},3),betaSNR(ia),'b','filled','MarkerFaceAlpha',0.5)

plot3([-70 -70+ofs],[60 60],[0 0],'k','LineWidth',2); hold on
plot3([-70 -70],[60 60+ofs],[0 0],'k','LineWidth',2)
plot3([-70 -70],[60 60],[0 0+ofs],'k','LineWidth',2)

view([-125 30]); 
xlim([-50 100]); ylim([-30 50]); zlim([0 80])
axis equal; box off; axis off

% Set 3
subplot(1,3,3)
% scatter3(locs(:,1),locs(:,2),locs(:,3),'k');
hold on
scatter3(roi(2,1),roi(2,2),roi(2,3),100,'r','filled')

scatter3(locs(indexSet{3},1),locs(indexSet{3},2),locs(indexSet{3},3),powDsync(ia),'b')
[c,ib] = intersect(indexSet{2},indexSet{4});
scatter3(locs(indexSet{4},1),locs(indexSet{4},2),locs(indexSet{4},3),powDsync(ib),'b','filled','MarkerFaceAlpha',0.5)

plot3([-70 -70+ofs],[60 60],[0 0],'k','LineWidth',2); hold on
plot3([-70 -70],[60 60+ofs],[0 0],'k','LineWidth',2)
plot3([-70 -70],[60 60],[0 0+ofs],'k','LineWidth',2)

view([-125 30]); 
xlim([-50 100]); ylim([-30 50]); zlim([0 80])
axis equal; box off; axis off

set(gcf,'Position',[680         585        1075         393])