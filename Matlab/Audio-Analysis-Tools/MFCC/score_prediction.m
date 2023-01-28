function [overall_accuracy, per_class_accuracy] = score_prediction(test_labels,predicted_labels)

%   Name: Karan Pareek
%   Student ID: kp2218
%
%
%   Compute the confusion matrix given the test labels and predicted labels.
%
%   Parameters
%   ----------
%   test labels: 1 x NE array
%       array of ground truth labels for test data predicted labels: 1 x NE test array
%       array of predicted labels
%
%   Returns
%   -------
%   overall accuracy: scalar
%         The fraction of correctly classified examples.
%   per class accuracy: 1 x 2 array
%       The fraction of correctly classified examples
%       for each instrument class.
%       per class accuracy[1] should give the value for
%       instrument class 1, per class accuracy[2] for
%       instrument class 2, etc.

%% Overall Accuracy

% Comparing the labels of both classes
total_accuracy = (predicted_labels == test_labels);

% Deriving the length of the correct labels
value_A = find(total_accuracy == 1);
value_A = length(value_A);

% Accuracy = Correct/Total
overall_accuracy = value_A/length(total_accuracy);

%% Per Class Accuracy

% Class 1
% Comparing the labels
l1 = find(test_labels == 1);
class1_accuracy = (test_labels(1:length(l1)) == predicted_labels(1:length(l1)));

% Deriving the length of the correct labels
value_1 = find(class1_accuracy == 1);
value_1 = length(value_1);

% Accuracy = Correct/Total
C1 = value_1/length(class1_accuracy);

% Class 2
% Comparing the labels
l2 = find(test_labels == 2);
class2_accuracy = (test_labels(1:length(l2)) == predicted_labels(1:length(l2)));

% Deriving the length of the correct labels
value_2 = find(class2_accuracy == 1);
value_2 = length(value_2);

% Accuracy = Correct/Total
C2 = value_2/length(class2_accuracy);

per_class_accuracy = [C1,C2];

end
