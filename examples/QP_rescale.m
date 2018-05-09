%% This script stores a solution of a section

thetas = atan(cpy(cpx<0.0001 & cpx > -0.0001 & cpy >= 0 & bdy ~=1)./(cpz(cpx<0.0001 & cpx > -0.0001 & cpy >= 0 & bdy ~=1)-cenf(3)));
sln = pert(cpx<0.0001 & cpx > -0.0001 & cpy >= 0 & bdy ~=1);

thetasdiag = atan(sqrt(cpx(abs(cpx - cpy) <0.0001 & cpy >= 0 & bdy ~=1).^2 + cpy(abs(cpx - cpy) <0.0001 & cpy >= 0 & bdy ~=1).^2)./(cpz(abs(cpx - cpy) <0.0001 & cpy >= 0 & bdy ~=1)-cenf(3)));
slndiag = pert(abs(cpx - cpy) <0.0001 & cpy >= 0 & bdy ~=1);

filename = 'QP_dx_05';

csvwrite(['data/', filename, '/theta.dat'],thetas);
csvwrite(['data/', filename, '/soln.dat'],sln);

csvwrite(['data/', filename, '/thetadiag.dat'],thetas);
csvwrite(['data/', filename, '/solndiag.dat'],sln);

fileID = fopen(['data/', filename, '/sln2.dat'],'wt');

fwrite(fileID, '[');

for i = 1:size(thetas,1)
    fwrite(fileID, ['[', num2str(thetas(i)),',',num2str(sln(i)),']']);
    if i ~= size(thetas,1)
       fwrite(fileID, [',', double(sprintf('\n'))]);
    end
end;

fwrite(fileID, ']');
fclose(fileID);

fileID2 = fopen(['data/', filename, '/sln2diag.dat'],'wt');

fwrite(fileID2, '[');

for i = 1:size(thetasdiag,1)
    fwrite(fileID2, ['[', num2str(thetasdiag(i)),',',num2str(slndiag(i)),']']);
    if i ~= size(thetasdiag,1)
       fwrite(fileID2, [',', double(sprintf('\n'))]);
    end
end;

fwrite(fileID2, ']');
fclose(fileID2);