function [eventTstart, eventTend] = mergeevents_time (eventTstart, eventTend, mergeTimeMax)
%[eventTstart, eventTend] = MERGEEVENTS_TIME ...
%                        (eventTstart, eventTend, mergeTimeMax)
%   Given times of event boudaries, merge events if the time interval
%   between them is inferior to mergeTimeMax
%
% INPUTS:
%   - eventTstart       : Starting time of each event (s)
%   - eventTend         : Ending time of each event (s)
%   - mergeTimeMax      : Maximum time duration for merging (s)
%
% OUTPUTS:
%   - eventTstart       : Merged starting time (s)
%   - eventTend         : Merged ending time (s)


nEvents         = length(eventTstart);
% nEventsOri      = nEvents;
i               = 1;
while i < nEvents-1
    % If events are close one from another (start of event i+1 close from
    % end of event i)
    if (eventTstart(i+1) - eventTend(i)) < mergeTimeMax
        % Merge the events (end of event i is now equal to end of event i+1)
        eventTend (i)      = eventTend (i+1);
        % Delete event i+1
        eventTstart (i+1)  = [];
        eventTend (i+1)    = [];  
        % Decrement nb_events
        nEvents = nEvents - 1;
    else
        i=i+1;
    end
end

end

