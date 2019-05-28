root = '\\nas.ls.berkeley.edu\111lab\Student-Redirect$\matthewnchow\My Documents\NVC\Data and Pictures\';
ddir = [root, 'testing_dir\'];
filebase = [ddir, 'Rabi_-5dBm_300ns_run'];

figure(1)
[dat, badruns] = PlotSimpRabi(filebase, 1:400);
length(badruns)
plot(dat(:,1), dat(:,2))

function [dat, badruns] = PlotSimpRabi(filebase, runs)
    d1 = importdata([filebase, num2str(runs(1)), '.csv']);
    ts = d1.data(:,1);
    cts = zeros(length(ts),1);
    dat = zeros(length(ts),2);
    badruns = [];
    for run = runs
       d = importdata([filebase, num2str(run), '.csv']);
       src = d.data(:,2);
       ref = d.data(:,3);
       if ~(sum(src) == 0 || sum(ref) == 0)
           cts = cts + src./ref;
       else 
           badruns = [badruns, run];
       end
    end
    dat(:,1) = ts;
    dat(:,2) = cts;
end