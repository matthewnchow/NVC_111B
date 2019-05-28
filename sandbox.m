% % Find a serial port object.
% MCLpiezo.fid = instrfind('Type', 'serial', 'Port', PortMap('PI_Piezo'), 'Tag', '');
% 
% % Create the serial port object if it does not exist
% % otherwise use the object that was found.
% if isempty(MCLpiezo.fid)
%     MCLpiezo.fid = serial(PortMap('PI_Piezo'),'BaudRate',115200,'DataBits',8);
% else
%     fclose(MCLpiezo.fid);
%     MCLpiezo.fid = MCLpiezo.fid(1);
% end

ddir = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\';
dfilename = 'test-data_5-18-19.csv';
d = importdata([ddir, dfilename]);


N = int32(500);
imax = idivide(length(d.data(:,1)),N)
xs = zeros(N,1);
ys = zeros(N,1);

for j = 1:N
    xs(j) = d.data(j,1);
    for i = 1:imax
        ys(j) = ys(j) + (d.data(j + (i-1)*N, 2)/double(imax));
    end
end

figure
plot(d.data(:,1), d.data(:,2), 'o')

figure 
plot(xs, ys)