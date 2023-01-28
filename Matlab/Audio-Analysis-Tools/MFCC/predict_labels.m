function predicted_labels = predict_labels(train_features,train_labels,test_features)

%   Name: Karan Pareek
%   Student ID: kp2218
%
%   Predict the labels of the test features,
%   given training features and labels,
%   using a nearest-neighbor classifier.
%
%   Parameters
%   ----------
%   train features: NF x NE train matrix
%       matrix of training set features (NF is number of
%       features and NE train is number of training feature instances) 
%   train labels: 1 x NE train array
%       vector of labels (class numbers) for each instance
%       of train features
%   test features: NF x NE test matrix
%       matrix of test set features (NF is number of
%       features and NE test is number of testing feature instances)
%
%   Returns
%   -------
%   predicted labels: 1 x NE test array
%       array of predicted labels

% Blank matrix for the predicted labels
predicted_labels = zeros(1,size(test_features,2));

% Computing the nearest neighbour algorithm
for n = 1:size(test_features,2)
    
    % Dot product of each test feature vector with all the feature vectors 
    % of the training data
    dot_product = train_features .* test_features(:,n);
    
    % Finding the maximum value of the derived vector
    m = max(dot_product);
    
    % Calculating the corresponding index(label) number
    [~,index] = max(m);

    % Assigning a label to the predicted label array
    predicted_labels(:,n) =+ train_labels(index);
    
end

end
