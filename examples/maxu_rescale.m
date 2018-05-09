%% This script rescales the maxu table to about 50 entries

filename = 'new_4915_dx_05_eps_2em6';

maxu = csvread(['data/', filename, '/maxu.dat']);
maxuscal = maxu(1:32:end,:);

csvwrite(['data/', filename, '/maxuscal.dat'],maxuscal);

fileID = fopen(['data/', filename, '/maxuscal2.dat'],'wt');

fwrite(fileID, '[');

for i = 1:size(maxuscal,1)
    fwrite(fileID, ['[', num2str(maxuscal(i,1)),',',num2str(maxuscal(i,2)),']']);
    if i ~= size(maxuscal,1)
       fwrite(fileID, [',', double(sprintf('\n'))]);
    end
end;

fwrite(fileID, ']');
fclose(fileID);