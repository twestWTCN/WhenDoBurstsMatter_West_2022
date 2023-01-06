function ABCErrorFX(R)
warning('Fitting failed')
try
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ERLIST'])
catch
    % errorlist didnt exist
    ERLIST = [];
end
LE = lasterror;
ERLIST{end+1,1} = [R.out.dag ' ' date ' ' LE.message];
save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ERLIST'],'ERLIST')