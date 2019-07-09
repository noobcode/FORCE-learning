seed_1 = 1; % seed of first run
seed_2 = 2; % seed of second run

load_weights_1 = 0; % initialize weights randomly in the first run
load_weights_2 = 1; % load weights in the second run

training_setting_1 = 0; % Nicola & Clopath in the first run
training_setting_2 = 3; % only evaluate in the second run

%% Nicola & Clopath training, load weights, replay
replay = 0; % ignore
IZFORCESINE(seed_1, load_weights_1, training_setting_1, replay);
IZFORCESINE(seed_2, load_weights_2, training_setting_2, replay);

%% Our approach
replay_1 = 1; % first run replay experiment
replay_2 = 2; % second run replay experiment
IZFORCESINE(seed_1, load_weights_1, training_setting_1, replay_1);
IZFORCESINE(seed_2, load_weights_2, training_setting_2, replay_2);