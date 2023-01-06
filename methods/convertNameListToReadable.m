function newLabList = convertNameListToReadable(R,labList)
for i = 1:numel(labList)
    listEl = labList(i);
    splitListEl = split(listEl,' ');
    
    tagWithoutNumeric = splitListEl{3}(1:regexp(splitListEl{3},'\d*','start')-1);
                        TagId = regexp(splitListEl{3},'\d*','Match');
TagId = str2num(TagId{1});
    switch splitListEl{1}
        case 'int'
            switch tagWithoutNumeric
                case 'T'
% %                     newLabList{i} = [splitListEl{2} ' time constant ' R.chsim_name{str2num(splitListEl{3}(end))}];
                    newLabList{i} = ['time constant ' R.chsim_name{str2num(splitListEl{3}(end))}];
                case 'alpha'
                    newLabList{i} = ['noise alpha ' R.chsim_name{str2num(splitListEl{3}(end))}];
                case 'Bn'
                    newLabList{i} = ['baseline firing ' R.chsim_name{str2num(splitListEl{3}(end))}];
                case 'Sn'
                    newLabList{i} = ['sig. slope ' R.chsim_name{str2num(splitListEl{3}(end))}];
                
                case 'C'
                    newLabList{i} = [R.chsim_name{str2num(splitListEl{3}(end))} ' input gain'];
                case 'G'
                    id = TagId;
                    if id == 1; newLabList{i} = ['mp -> mp']; end
                    if id == 2; newLabList{i} = ['mp -> sp']; end
                    if id == 3; newLabList{i} = ['ii -> mp']; end
                    if id == 4; newLabList{i} = ['ii -> ii']; end
                    if id == 5; newLabList{i} = ['mp -> ii']; end
                    if id == 6; newLabList{i} = ['dp -> ii']; end
                    if id == 7; newLabList{i} = ['sp -> sp']; end
                    if id == 8; newLabList{i} = ['sp -> mp']; end
                    if id == 9; newLabList{i} = ['ii -> dp']; end
                    if id == 10; newLabList{i} = ['dp -> dp']; end
                    if id == 11; newLabList{i} = ['sp -> dp']; end
                    if id == 12; newLabList{i} = ['ii -> sp']; end
                    if id == 13; newLabList{i} = ['sp -> ii']; end
                    if id == 14; newLabList{i} = ['dp -> sp']; end
                    
                otherwise
                    warning('List element not recognized')
                    
            end
    end
    
end
