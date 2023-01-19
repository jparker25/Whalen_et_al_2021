figure
subplot(1,2,1)
hold on
%plot([0,-cent_int/cent_slope],[cent_int,0],'--','Color',[.9 .9 .9],'LineWidth',5)
scatter(twoHz.cents(1,1), twoHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(twoHz.cents(2,1), twoHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(twoHz.non_osc_cent(1), twoHz.non_osc_cent(2), 2000, 'k', '.')

scatter(twoHz.xnoG,twoHz.ynoS,200,[0.8 0.8 0.8],'.')
scatter(twoHz.strG_total(twoHz.pttype_osc==-1),twoHz.strS_total(twoHz.pttype_osc==-1),200,'b','.')
scatter(twoHz.strG_total(twoHz.pttype_osc==1),twoHz.strS_total(twoHz.pttype_osc==1),200,'r','.')

plot([0 800],[400 100],'k','LineWidth',3)
line([400 400],[0 250],'Color','k','LineWidth',3)
hold off

[hLg, icons]=legend({'AP','IP',"Non-Osc"})
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons(1:3),'MarkerSize',20);
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])

subplot(1,2,2)
hold on
%plot([0,-cent_int/cent_slope],[cent_int,0],'--','Color',[.9 .9 .9],'LineWidth',5)
scatter(fifteenHz.cents(1,1), fifteenHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(fifteenHz.cents(2,1), fifteenHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(fifteenHz.non_osc_cent(1), fifteenHz.non_osc_cent(2), 2000, 'k', '.')
scatter(fifteenHz.xnoG,fifteenHz.ynoS,200,[0.8 0.8 0.8],'.')
plot([0 800],[400 100],'k','LineWidth',3)
line([400 400],[0 250],'Color','k','LineWidth',3)
for type = [1 -1]
    col = cols(1+(type<0),:);
    scatter(fifteenHz.strG_total(fifteenHz.pttype_osc==type),fifteenHz.strS_total(fifteenHz.pttype_osc==type),200,col,'.')
end

hold off
[hLg, icons]=legend({'AP','IP',"Non-Osc"})
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons(1:3),'MarkerSize',20);
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])


figure;
hold on
fifteen_all = []
for i = 1:length(fifteenHz.xnoG)
    fifteen_all = [fifteen_all; [fifteenHz.xnoG(i) fifteenHz.ynoS(i);]];
end

xx15 = fifteenHz.strG_total(fifteenHz.pttype_osc==-1);
yy15 = fifteenHz.strS_total(fifteenHz.pttype_osc==-1);
for i = 1:length(xx15)
    fifteen_all = [fifteen_all; [xx15(i) yy15(i);]];
end

xx15 = fifteenHz.strG_total(fifteenHz.pttype_osc==1);
yy15 = fifteenHz.strS_total(fifteenHz.pttype_osc==1);
for i = 1:length(xx15)
    fifteen_all = [fifteen_all; [xx15(i) yy15(i);]];
end

scatter(twoHz.cents(1,1), twoHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(twoHz.cents(2,1), twoHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(mean(fifteen_all(:,1)),mean(fifteen_all(:,2)),2000,[0.9 0.9 0.9],'.')
scatter(fifteen_all(:,1),fifteen_all(:,2),200,[0.7 0.7 0.7],'.')

%scatter(twoHz.non_osc_cent(1), twoHz.non_osc_cent(2), 2000, 'k', '.')
%scatter(twoHz.xnoG,twoHz.ynoS,200,[0.8 0.8 0.8],'.')

scatter(twoHz.strG_total(twoHz.pttype_osc==-1),twoHz.strS_total(twoHz.pttype_osc==-1),200,'b','.')
scatter(twoHz.strG_total(twoHz.pttype_osc==1),twoHz.strS_total(twoHz.pttype_osc==1),200,'r','.')



%scatter(fifteenHz.cents(1,1), fifteenHz.cents(1,2), 2000, [0.9 0.9 0.9], 'x') 
%scatter(fifteenHz.cents(2,1), fifteenHz.cents(2,2), 2000, 'k', 'x')
%scatter(fifteenHz.non_osc_cent(1), fifteenHz.non_osc_cent(2), 2000, 'k', 'x')



%scatter(fifteenHz.xnoG,fifteenHz.ynoS,200,[0.8 0.8 0.8],'x')
%scatter(fifteenHz.strG_total(fifteenHz.pttype_osc==-1),fifteenHz.strS_total(fifteenHz.pttype_osc==-1),200,'k','.')
%scatter(fifteenHz.strG_total(fifteenHz.pttype_osc==1),fifteenHz.strS_total(fifteenHz.pttype_osc==1),200,[0.9 0.9 0.9],'.')

hold off

legend({'AP','IP',"15Hz"})
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])

