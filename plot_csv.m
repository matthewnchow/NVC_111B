f = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\150um_5-17_FSpectra.csv';

d = importdata(f);
x = d.data(:,1);
y = d.data(:,2);
figure
hax = axes;
plot(x,y);
xlabel(strcat(d.colheaders(1), ' (nm)'));
ylabel(d.colheaders(2));
line([637 637],get(hax,'YLim'),'Color',[1 0 0])
line([575 575],get(hax,'YLim'),'Color',[1 1 0])