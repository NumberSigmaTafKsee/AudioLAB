function [train_features, train_labels] = create_set(fpath1, fpath2, params)

%   Name: Karan Pareek
%   Student ID: kp2218
%
%   Compute features and parameters for training data.
% 
%   Parameters
%   ----------
%   fpath1: string
%       full path to audio file with training data from class 1
%   fpath2: string
%       full path to audio file with training data from class 2
%   params: struct
%       Matlab structure with fields are win size, hop size,min freq, max freq, 
%       num mel filts, n dct, the parameters needed for computation of MFCCs
%   Returns
%   -------
%   features: NF x NE matrix
%       matrix of training/testing set features (NF is number of
%       features and NE is number of feature instances)
%   labels: 1 x NE array
%       vector of training/testing labels (class numbers) for each instance
%       of features


[mfccs_1] = compute_mfccs(fpath1,params.win_size,params.hop_size,...
            params.min_freq,params.max_freq,params.num_mel_filts,params.n_dct);
[mfccs_2] = compute_mfccs(fpath2,params.win_size,params.hop_size,...
            params.min_freq,params.max_freq,params.num_mel_filts,params.n_dct);

% NF - number of features (columns)
% NE - number of feature instances (rows)

% Feature matrix (NF x NE)
train_features = [mfccs_1,mfccs_2];

% Traning/testing labels (1 x NE)
train_labels = [ones(1,size(mfccs_1,2)),(2*ones(1,size(mfccs_2,2)))];

% labels = [1,1,1,1,...,1,2,2,2,2,...,2]

end
