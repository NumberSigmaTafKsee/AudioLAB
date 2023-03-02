%{
  MP2
  GROUP 8:
    - CRUZ, Airon John
    - HERNANDEZ, Pierre Vincent
    ASSIGNED EQUATION:
                { +k,   when -pi   < x < -pi/2
    8.) f(x) =  {  0,   when -pi/2 < x <  pi/2
                { +k,   when  pi/2 < x <  pi
%}



% ====== START OF PROGRAM ======

% Display periodic function to be transformed to its Fourier Series
fprintf("\nThis script plots the periodic function f(x) and its Fourier series partial sum.\n");
fprintf("f(x) = +k,   when -pi   < x < -pi/2\n");
fprintf("f(x) =  0,   when -pi/2 < x <  pi/2\n");
fprintf("f(x) = +k,   when  pi/2 < x <  pi\n\n");


% ====== END OF PROGRAM  =======


% ====== FUNCTIONS =============

function [user_input] = get_num_input(strPrompt)
%{
    This function handles user's numerical input and returns
    it when it is a valid input.
    Parameter:
        strPrompt --> the String indicates the instruction of the program
                        for the respective asked values from the users.
%}
    tmp = '';

    % Continue asking user if input is empty
    % or not a computable input
    while(isempty(tmp))
        tmp = input(strPrompt);

        % Prompt user for empty input
        if(isempty(tmp))
            fprintf("Empty input is not allowed!\n\n");
        end
    end

    %this assigns the value of tmp to the return value of user_input
    user_input = tmp;
end

function [k] = get_amplitude()
%{
    This function calls for the get_num_input function to accept input
    values for preferred amplitude of the user. In addition, this function
    will add specific restrictions in order to handle such instances
    wherein desired inputs are not applicable.
    Parameters: NONE
%}
    k = 0;
    while k < 1
       k = get_num_input('Enter the amplitude parameter, k: ');

       % Prompt user for invalid input
       if(k < 1)
           fprintf("Input must be greater than 0!\n\n");
       end
    end
end

function [count] = get_sinusoid_count()
%{
    This function will used get_num_input function to ask for user's
    preferred sinusoid count input. This function will also implement
    restrictions specific to this function such that the program will
    continuously asks user for a valid input value.

    Parameters: NONE
%}
    count = 0;
    while count < 1
       count = get_num_input('Enter number of sinusoids to add: ');

       % Prompt user for invalid input
       if(count < 1)
           fprintf("Input must be greater than 0!\n\n");
       end
    end
end

function [fourier] = fourier_series(amplitude, sinusoidCount, xValues)
%{
    This function gets the fourier solution for each values of
    `x` (xValues) with the specified `amplitude` and `sinusoidCount`
    (number of sinusoids in the summation).
    Parameters:
        amplitude     --> the desired k-value of the user
        sinusoidCount --> the user's desired number of sinusoid
        xValues       --> the generated X values (Abscissa).
%}
    % Maximum number of elements for fourier.
    % Same number of elements as xValues
    maxi = numel(xValues);

    % Allocate memory for fourier with the same size as xValues
    fourier = zeros(size(xValues));

    % Get the constant
    a_0 = get_constant(amplitude);

    % Get the coefficients
    a_n = get_coefficients(amplitude, sinusoidCount);

    % Create symoblic scalar variables `x` and `n`
    syms x n
    % Create the fourier series function with the specified
    % number of sinusoids
    f(x) = a_0 + sum(a_n.*subs(cos(n*x), n, 1:sinusoidCount));

    % Loop through every element of `xValues`
    for i = 1:maxi
        % Store the answer of the `f(x)` for each element of `xValues`
        fourier(i) = f(xValues(i));
    end
end

function [a_0] = get_constant(amplitude)
%{
    This function will compute for the a_0 of the fourier series
    (constant) given an amplitude from the user's input.

    Parameters:
        amplitude --> the desired k-value of the user
%}

% this is derived from the analytical solution that is manually performed
    a_0 = amplitude/2;
end

function [a_n] = get_coefficients(amplitude, sinusoidCount)
%{
    This function will compute for the coefficients (a_n)
    of the sinusoids of the fourier series
    Parameters:
        amplitude     --> the user input for the desired amplitude.
        sinusoidCount --> the user's desired number of sinusoid.
%}
    % Generate n-values from 1 to sinusoidCount
    n = 1:1:sinusoidCount;

    % Compute for the coefficients
    a_n = (-2*amplitude)./(pi.*n).*sin(n.*pi./2);
end

function plot_solution(amplitude, sinusoidCount)
%{
    This function plots the solution of the fourier series
    on top of its equivalent periodic function.
   Parameters:
        amplitude     --> the user input for the desired amplitude.
        sinusoidCount --> the user's desired number of sinusoid.
%}
    % Generate x-values
    xVals = linspace(-pi,pi,200);

    % Values for periodic function
    fx = [amplitude*ones(1,50),0*ones(1,100),amplitude*ones(1,50)];

    % Values/Solution for fourier series
    fourier = fourier_series(amplitude, sinusoidCount, xVals);

    % Setup figure
    figure, clf; hold on

    % Plot the periodic function, f(x)
    plot(xVals,fx);

    % Superimpose plot of Fourier Series
    plot(xVals,fourier);

    % Add title and labels
    xlabel('x');
    ylabel('f(x)');
    title('Plot of f(x) and its Fourier Series');
    hold off
end

% Get amplitude from user
amplitude = get_amplitude();

% Get number of sinusoids from user
sinusoidCount = get_sinusoid_count();

% Output/Display plot
plot_solution(amplitude, sinusoidCount);

