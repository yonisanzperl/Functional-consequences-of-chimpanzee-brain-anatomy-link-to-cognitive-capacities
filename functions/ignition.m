function stdevokedintegS=ignition(events,N,Tmax)

Isubdiag = find(tril(ones(N),-1));
nTRs = 5; %
T = 1:Tmax; 

%%% integration
% obtain 'events connectivity matrix' and integration value (integ)
% for each time point
for t = T
    for i = 1:N
        for j = 1:N
            phasematrix(i,j) = events(i,t)*events(j,t);
        end
    end
    cc = phasematrix; %*Cbin;
    cc = cc-eye(N);

    [comps csize] = get_components(cc);
    integ(t) = max(csize)/N;
end % end obtain integ

%%%% event trigger
nevents2 = zeros(1,N);
% save events and integration values for nTRs after the event
for seed = 1:N
    flag = 0;
    for t = T
        % detect first event (nevents = matrix with 1xnode and number of events in each cell)
        if events(seed,t) == 1 && flag == 0  % if events(seed,t-9) == 1 && flag == 0
            flag = 1;
            % events for each subject
            nevents2(seed) = nevents2(seed)+1;
        end
        % save integration value for nTRs after the first event (nodesx(nTR-1)xevents)
        if flag > 0
            % integration for each subject
            IntegStim2(seed,flag,nevents2(seed)) = integ(t);
            flag = flag+1;
        end
        % after nTRs, set flag to 0 and wait for the next event (then, integ saved for nTRs -1 events)
        if flag == nTRs
            flag = 0;
        end
    end
end

% std of the max ignition in the nTRs for each subject and for each node
for seed = 1:N
    stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
end

% std ignition across events for each subject in each node(Single Subject, S)
stdevokedintegS(:) = stdevokedinteg2;