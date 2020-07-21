clc; clear;

load('grainsize.mat')
numSieve = 10;
scale = .10769;
sieveSz = linspace(0,45/scale,numSieve); %size of sieve in equal intervals
files = fieldnames(repos);
pass = zeros(numSieve,length(files));

totVol = 0;
for sieve = 1:numSieve
    for agg = 1:length(files)
        grainsize = repos.(files{agg,1}).GrainSize;
        volume = repos.(files{agg,1}).Volume;
        if grainsize < sieveSz(sieve)
            pass(sieve,agg) = volume;
        end
        
        if sieve == 1
            totVol = totVol + volume;
        end
    end
end

figure(1)
results = (sum(pass,2)./totVol) * 100;
semilogx(sieveSz,results,'-*')
ylim([0 100])
xlim([0 45/scale])
title('Grain Size Distribution')
xlabel('Grain Size (pixels)')
ylabel('Cummulative Passing (%), Volume')
grid on