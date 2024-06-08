% Base path to your BIDS dataset
dataPath = '/home/mcesped/scratch/Datasets/bids_ieeg/Subj_space/';

% Get a list of all subjects
subjects = dir(fullfile(dataPath, 'sub-*'));
subjects = subjects([subjects.isdir]);

% Loop over each subject
for i = 1:length(subjects)
    subjectPath = fullfile(dataPath, subjects(i).name);

    % Get a list of all sessions for the subject
    sessions = dir(fullfile(subjectPath, 'ses-*'));
    sessions = sessions([sessions.isdir]);

    % Loop over each session
    for j = 1:length(sessions)
        sessionPath = fullfile(subjectPath, sessions(j).name, 'ieeg');

        % Get a list of EEG files
        eegFiles = dir(fullfile(sessionPath, '*.edf')); 

        % Loop over each EEG file
        for k = 1:length(eegFiles)
            eegFilePath = fullfile(sessionPath, eegFiles(k).name);

            % Load the EEG data
            cfg = [];
            cfg.dataset = eegFilePath;
            eegData = ft_preprocessing(cfg);

            % Save the processed data
            save(fullfile(sessionPath, ['processed_' eegFiles(k).name]), 'eegData');
        end
    end
end
