function [montage, bip_data] = getbipolarmontage( labels, special, data )
%Create bipolar for HFO and spike detection
Channel_labels_Raw = labels;
% Find element which is not number in each label
% Channel_labels_Init = cellfun(@(x) cell2mat(regexp(x,'[^0-9]','match')),Channel_labels_Raw);
Channel_labels_Init = cellfun(@(x) cell2mat(regexp(x,'[^0-9]','match')),Channel_labels_Raw,'UniformOutput',0);
empty_set = find(cellfun(@isempty, Channel_labels_Init));
if ~isempty(empty_set),
    Channel_labels_Init = Channel_labels_Init(~cellfun(@isempty, Channel_labels_Init));
end
Electrode_num = length(unique(Channel_labels_Init));
[Electrode_Name,~,b] = unique(Channel_labels_Init,'stable');
% Initialize the montage matrix
tra = zeros(length(labels) - Electrode_num,length(labels));
labelorg = Channel_labels_Raw;
labelnew = cell(length(labelorg)-Electrode_num,1)';
% create montage matrix
Channel_ind_new = 1;
for i = 1:Electrode_num
%     Electrode_Contacts = length(b(b == i));
    tmp = regexp(labelorg(b==i), '[0-9]{1,}', 'match');
    Electrode_Contacts =  max(cellfun(@str2double, tmp));
    if Electrode_Contacts > 99,
        str = '%03i';
    else
        str = '%02i';
    end
    for j = 1:Electrode_Contacts
        if exist('special','var')&~isempty(special)&strcmp(special{1}, Electrode_Name{i})
            % special indexing
            if mod(j, special{2})==0,
                Channel_ind_new = Channel_ind_new + 1;
                continue;
            end
        elseif Electrode_Contacts > 30 & mod(j, 8)==0, % must be a grid
                % if grid, then don't take differences between 8-9, etc..
                Channel_ind_new = Channel_ind_new + 1;
                continue;
        end
        name1 = strcat(Electrode_Name{i},num2str(j,str));
        name2 = strcat(Electrode_Name{i},num2str(j+1,str));
        % locate this two channels
        [~,Locb1] = ismember(name1,Channel_labels_Raw);
        [~,Locb2] = ismember(name2,Channel_labels_Raw);
        
        
        
        if Locb1 > 0 & Locb2 > 0,
            labelnew{1,Channel_ind_new}=sprintf('%s-%s',name1,name2);
            tra(Channel_ind_new,Locb1) = 1;
            tra(Channel_ind_new,Locb2) = -1;
        end
        Channel_ind_new = Channel_ind_new + 1;
    end
end
emp = find(cellfun(@isempty, labelnew));
tra(emp, :) = [];
labelnew(emp) = [];
montage.tra = tra;
montage.labelnew = labelnew;
montage.labelorg = labelorg;

if nargout>1,
    % create the new dataset
    bip_data = [];
    for j=1:size(montage.tra,1),
        ch1 = find(montage.tra(j, :)==1);
        ch2 = find(montage.tra(j, :)==-1);
        bip_data(:, j) = data(:, ch1)-data(:, ch2);
    end

end