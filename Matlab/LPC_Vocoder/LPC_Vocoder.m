classdef LPC_Vocoder < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        LoadSoundFile_Button   matlab.ui.control.Button
        OrigSignal_Figure      matlab.ui.control.UIAxes
        Title_Label            matlab.ui.control.Label
        EncodeDecode_Button    matlab.ui.control.Button
        SynthSignal_Figure     matlab.ui.control.UIAxes
        Author_Label           matlab.ui.control.Label
        BounLogo_Image         matlab.ui.control.Image
        PlayOrigSound_Button   matlab.ui.control.Button
        PlaySynthSound_Button  matlab.ui.control.Button
        Save_Button            matlab.ui.control.Button
    end

    
    properties (Access = private)
        orig_signal; %A variable to hold the original input speech signal
        synth_signal; %A variable to hold the synthesized output speech signal
        fs; %A variable to hold the original sampling frequency
        order = 10; %A variable to hold the prediciton order of the LPC algorithm 
    end
    
    methods (Access = private)
        
        function [coeffs, pitches, voices, gains] = LPC_Encoder(app, orig_signal, fs, M)
            %This function encodes the given original speech signal (with sampling frequency fs) with Linear Predictive 
            %Coding Model as specified in page 269 of "Speech Coding Algorithms" textbook of Wai Chu. 
            
            
            %Returns the parameters of the model:
            %   Excitation parameters:
            %       -isVoiced = Voiced/Unvoiced Classification for each frame
            %       -pitches = Pitch periods for each voiced speech frame 
            %       -gains = Gains for each speech frame
            %   Vocal tract system parameters
            %       -coeffs = Coefficients (ai, i = 1,2,...,M) of the all-pole digital filter for each frame 
            
            %1)Initializing necessary variables to be used during algorithm.
            frame_size = 30e-3; %In units of seconds
            frame_length = round(fs .* frame_size); %In units of samples
            
            %2)Get voiced/unvoiced classification and pitch period of frames.
            [voices, pitches] = voicing_detector(app, orig_signal, fs, frame_length);
            
            %3) Frame segmentation and calc of coeffs and gains
            for index = 1 : frame_length : (length(orig_signal) - frame_length)
                frame_signal = filter([1 -0.9378], 1, orig_signal(index : (index + frame_length -1))); %Obtain one frame of samples and store it in an array as the current frame.
                %Coeffs:
                [frame_coeffs, frame_coeffs_count, error] = levinson_durbin(app, frame_signal, M); %Calculates the coeffs for the current frame using the Levinson-Durbin algorithm.
                coeffs(index : (index + frame_coeffs_count - 1)) = frame_coeffs; %Store the frame coeffs inside the array of all coeffs.
                %Gain:
                frame_pitch = pitches(index);
                frame_voice = voices(index);
                gains(index) = gain_calc(app, error, frame_voice, frame_pitch);
            end
        end
        
        function synth_signal = LPC_Decoder(app, coeffs, pitches, voices, gains)
            %This function decodes the given Linear Predictive Coding Model
            %parameters and sytnhesizes a speech signal based on the model for 
            %speech production as specified in the page 475 of "Theory and
            %Applications of Digital Speech Processing" textbook by L. Rabiner and R. Schafer
            %and as specified in the page 271 of "Speech Coding Algorithms" textbook of Wai Chu.
            
            %Returns:
            %  synth-signal = the synthesized speech signal 
            
            frame_length=1;
            for i=2:length(gains)
                if gains(i) == 0
                    frame_length = frame_length + 1;
                else 
                    break;
                end
            end

            for index = 1 : frame_length : (length(gains))   
                if voices(index) == 1
                    frame_pitch = pitches(index);
                    synth_frame_signal = synth_voiced(app, coeffs, gains, frame_length, frame_pitch, index);
                else 
                    synth_frame_signal = synth_unvoiced(app, coeffs, gains, frame_length, index); 
                end
                synth_signal(index : index + frame_length - 1) = synth_frame_signal;
            end
        end
        
        function synth_frame_signal = synth_unvoiced(app, coeffs, gains, frame_length, index)
            %This function synthesizes unvioced speechs based on the speech production model as specified 
            %in the page 475 of "Theory and Applications of Digital Speech Processing" textbook by L. Rabiner and R. Schafer
            white_noise = randn(1, frame_length);
            synth_frame_signal = filter(1, [1 coeffs((index + 1) : (index + 1 + 9))], white_noise) .* gains(index); %9 comes from the prediction order as 10 - 1
        end
        
        function synth_frame_signal = synth_voiced(app, coeffs, gains, frame_length, frame_pitch, index)
            %This function synthesizes vioced speechs based on the speech production model as specified 
            %in the page 475 of "Theory and Applications of Digital Speech Processing" textbook by L. Rabiner and R. Schafer
            for f = 1 : frame_length
                if f ./ frame_pitch == floor(f ./ frame_pitch)
                    pulse_train(f) = 1;
                else 
                    pulse_train(f) = 0;
                end
            end
            synth_frame_signal = filter(1, [1 coeffs((index + 1) : (index + 1 + 9))], pulse_train) .* gains(index);
        end
        
        function [voices, pitches] = voicing_detector(app, orig_signal, fs, frame_length)
            %This funtion classifies the frames of the signal as voiced
            %or unvoiced based on the principles explained in page 271 of 
            %"Speech Coding Algorithms" textbook of Wai Chu.
            
            %Parameters:
            %   -orig_signal = sound signal
            %   -frame_length = is the length of frame in units of samples
            
            %Returns:
            %   -voices = an array pf 0s and 1s indicarting voicing of the frame
            %   -pitcghes = an array of pitch periods 
            
            for index = 1 : frame_length : (length(orig_signal) - frame_length) %Calculating magnitude sum energy, zero crossing rate and pitch period parameters for each frame
                frame_signal = filter([1 -0.9378], 1, orig_signal(index : index + frame_length - 1));
                
                ms_energy(index : (index + frame_length - 1)) = ms_energy_calc(app, frame_signal);
                zc_rate(index : (index + frame_length - 1)) = zc_rate_calc(app, frame_signal);
                pitches(index : (index + frame_length - 1)) = pitch_calc(app, frame_signal, fs);
            end
            
            %Obtaining thresholds for magnitude sum energy, zero crossing
            %rate and pitch period thresholds and classifying each frame based on these metrics 
            thresh_ms_energy = (((sum(ms_energy) ./ length(ms_energy)) - min(ms_energy)) .* 0.67) + min(ms_energy);
            voiced_ms_energies = ms_energy > thresh_ms_energy; 
            
            thresh_zc_rate = (((sum(zc_rate) ./ length(zc_rate)) - min(zc_rate)) .* 1.5) + min(zc_rate);
            voiced_zc_rates = zc_rate < thresh_zc_rate;
            
            thresh_pitch = (((sum(pitches) ./ length(pitches)) - min(pitches)) .* 0.5) + min(pitches);
            voiced_pitches = pitches > thresh_pitch;
            
            %Combining seperate classifications for each parameter and
            %deciding on whether the frame is voiced or unvoiced
            for index = 1 : (length(orig_signal) - frame_length)
                if voiced_ms_energies(index) .* voiced_zc_rates(index) .* voiced_pitches(index) == 1
                    voices(index) = 1;
                else
                    voices(index) = 0;
                end
            end
        end
        
        function ms_energy = ms_energy_calc(app, frame_signal)
            %This function calculates the magnitude sum energy of a signal
            %frame based on the equations on page 272 of "Speech Coding Algorithms" textbook of Wai Chu.
            
            [B, A] = butter(9, 0.33, "low");
            ms_energy = sum(abs(filter(B, A, frame_signal))); % Voiced speech has energy concentrated in low-freq region due to the relatively low value of pitch freq. 
                                                              % Better discrimination can be achieved through lowpass filtering.
        end
        
        function zc_rate = zc_rate_calc(app, frame_signal)
            %This function calculates the zero crossing rate of a signal
            %frame based on the equations on page 272 of "Speech Coding Algorithms" textbook of Wai Chu.
            
            zc_rate = 0;
            for n = 1 : length(frame_signal)
                if n + 1 > length(frame_signal)
                    break
                end
                zc_rate = zc_rate + (1 ./ 2) .* abs(sign(frame_signal(n+1)) - sign(frame_signal(n)));
            end      
        end
        
        function pitch = pitch_calc(app, frame_signal, fs)
            %This function calculates the pitch period of a signal
            %frame based on the equations on page 34 of "Speech Coding Algorithms" textbook of Wai Chu.
            min_period = round(fs .* 2e-3);
            max_period = round(fs .* 20e-3);
            
            R = xcorr(frame_signal);
            [R_max, R_mid] = max(R);
     
            pitch_per_range = R(R_mid + min_period : R_mid + max_period);
            [R_max, R_mid] = max(pitch_per_range);
     
            pitch = R_mid + min_period;      
        end
        
        function [frame_coeffs, frame_coeffs_count, error] = levinson_durbin(app, frame_signal, M)
            %This function calculates the vocal tract system parameters and
            %error for the Linear Predictive Coding Model for the given
            %frame signal based on the Levinson-Durbin Method (with prediction order M)
            %as specified in page 112 of "Speech Coding Algorithms" textbook of Wai Chu.
            
            %Arguments:
            %   -frame_signal = the given frame of speech signal
            %   -M = prediction order to be used in the method
            
            %Returns:
            %   -frame_coeffs = array of prediction coefficients of the all-pole digital filter.
            %   -frame_coeffs_count = size of the frame_coeffs array
            %   -error = difference between the original frame signal and
            %   the estimated frame signal based on the frame coeffs
            
            temp_coeffs = [zeros(M+1); zeros(M+1)]; %Initialize a matrix of zeros to record coefficients later, this matrix is updated all the time
            
            %Get autocorrelation coefficients matrix
            temp_auto_corr = xcorr(frame_signal); %Calculate to autocorrelation of the frame signal and record it in a temporary variable.
            R = temp_auto_corr(((length(temp_auto_corr) + 1) ./ 2) : length(temp_auto_corr)); %R is the array of R[l] 
            
            % Manually get predictor parameter for the 0 order
            %step = 1;
            J(1) = R(1);
            
            %Get predictor parameters for the higher orders
            for step = 2 : (M + 1)
                s_k = 0; %Summation term, set to 0 for each iteration
                for i = 2 : (step - 1)
                    s_k = s_k + temp_coeffs(i, (step - 1)) .* R(step - i + 1); 
                end
                
                k(step) = (R(step) + s_k) ./ J(step - 1);
                J(step) = J(step - 1) .* (1 - k(step) .^ 2);
                
                temp_coeffs(step, step) = -1 * k(step);
                temp_coeffs(1, step) = 1;
                
                for i = 2 : (step - 1)
                    temp_coeffs(i, step) = temp_coeffs(i, (step - 1)) - k(step) .* temp_coeffs((step - i + 1), (step - 1));
                end
            end
            
            frame_coeffs = temp_coeffs((1:step), step)';
            frame_coeffs_count = length(frame_coeffs);
            
            %Prediction error calculation
            est_frame_sig = filter([0 -frame_coeffs(2:end)], 1, frame_signal);
            error = frame_signal - est_frame_sig; 
        end
        
        function [gain, power] = gain_calc(app, error, frame_voice, frame_pitch)
            %This function calculates the gain of a given frame of speech signal based on the 
            %algorithm provided in the in page 270 of "Speech Coding Algorithms" textbook of Wai Chu.
            
            %Argumnents:
            %   -error = difference between the original signal and the approximated one
            %   =frame_voice = 0 if frame is unvioced, 1 if frame is voiced
            %   -frame_pitch = pitch of the frame
            
            %Returns:
            %   -gain = gain of signal in the frame
            %   -power = power of the signal in the frame
            
            if frame_voice == 0 %Unvoiced Frame
                denom = length(error);
                power = sum(error(1:denom) .^ 2) ./ denom;
                gain = sqrt(power);
            else %Voiced Frame
                denom = (floor(length(error) ./ frame_pitch) .* frame_pitch);
                power = sum(error(1:denom) .^ 2) ./ denom;
                gain = sqrt(power .* frame_pitch);
            end
        end
      
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadSoundFile_Button
        function LoadSoundFile_ButtonPushed(app, event)
            filename = uigetfile("MultiSelect","off");
            [app.orig_signal, app.fs] = audioread(filename); 
            plot(app.OrigSignal_Figure, app.orig_signal);
        end

        % Button pushed function: EncodeDecode_Button
        function EncodeDecode_ButtonPushed(app, event)
            [coeffs, pitches, voices, gains] = LPC_Encoder(app, app.orig_signal, app.fs, app.order); %ENCODE
            app.synth_signal = LPC_Decoder(app, coeffs, pitches, voices, gains); %DECODE
            plot(app.SynthSignal_Figure, app.synth_signal);
        end

        % Button pushed function: PlayOrigSound_Button
        function PlayOrigSound_ButtonPushed(app, event)
            soundsc(app.orig_signal, app.fs) 
        end

        % Button pushed function: PlaySynthSound_Button
        function PlaySynthSound_ButtonPushed(app, event)
            soundsc(app.synth_signal, app.fs)
        end

        % Button pushed function: Save_Button
        function Save_ButtonPushed(app, event)
           filename = inputdlg({'Enter a filename:'},'Save Synthesized Sound',[1 35],{'test_output_1.wav'});
           audiowrite(filename{1}, app.synth_signal, app.fs);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create LoadSoundFile_Button
            app.LoadSoundFile_Button = uibutton(app.UIFigure, 'push');
            app.LoadSoundFile_Button.ButtonPushedFcn = createCallbackFcn(app, @LoadSoundFile_ButtonPushed, true);
            app.LoadSoundFile_Button.BackgroundColor = [0 1 0];
            app.LoadSoundFile_Button.FontName = 'Comic Sans MS';
            app.LoadSoundFile_Button.FontSize = 15;
            app.LoadSoundFile_Button.FontWeight = 'bold';
            app.LoadSoundFile_Button.FontColor = [0.149 0.149 0.149];
            app.LoadSoundFile_Button.Position = [33 299 163 30];
            app.LoadSoundFile_Button.Text = 'LOAD SOUND FILE';

            % Create OrigSignal_Figure
            app.OrigSignal_Figure = uiaxes(app.UIFigure);
            title(app.OrigSignal_Figure, 'Original Speech Signal')
            xlabel(app.OrigSignal_Figure, 'Samples')
            ylabel(app.OrigSignal_Figure, 'Amplitude')
            app.OrigSignal_Figure.FontName = 'Comic Sans MS';
            app.OrigSignal_Figure.Position = [222 206 300 185];

            % Create Title_Label
            app.Title_Label = uilabel(app.UIFigure);
            app.Title_Label.HorizontalAlignment = 'center';
            app.Title_Label.FontName = 'Comic Sans MS';
            app.Title_Label.FontSize = 22;
            app.Title_Label.FontWeight = 'bold';
            app.Title_Label.Position = [112 409 436 34];
            app.Title_Label.Text = 'Linear Predictive Coding (LPC) Vocoder ';

            % Create EncodeDecode_Button
            app.EncodeDecode_Button = uibutton(app.UIFigure, 'push');
            app.EncodeDecode_Button.ButtonPushedFcn = createCallbackFcn(app, @EncodeDecode_ButtonPushed, true);
            app.EncodeDecode_Button.BackgroundColor = [0 1 0];
            app.EncodeDecode_Button.FontName = 'Comic Sans MS';
            app.EncodeDecode_Button.FontSize = 15;
            app.EncodeDecode_Button.FontWeight = 'bold';
            app.EncodeDecode_Button.FontColor = [0.149 0.149 0.149];
            app.EncodeDecode_Button.Position = [33 128 163 30];
            app.EncodeDecode_Button.Text = 'ENCODE & DECODE';

            % Create SynthSignal_Figure
            app.SynthSignal_Figure = uiaxes(app.UIFigure);
            title(app.SynthSignal_Figure, 'Synthesized Speech Signal with LPC')
            xlabel(app.SynthSignal_Figure, 'Samples')
            ylabel(app.SynthSignal_Figure, 'Amplitude')
            app.SynthSignal_Figure.FontName = 'Comic Sans MS';
            app.SynthSignal_Figure.Position = [222 13 300 185];

            % Create Author_Label
            app.Author_Label = uilabel(app.UIFigure);
            app.Author_Label.HorizontalAlignment = 'center';
            app.Author_Label.Position = [542 406 77 41];
            app.Author_Label.Text = {'Utku Türkbey'; 'BOUN EE'; 'V1.0'};

            % Create BounLogo_Image
            app.BounLogo_Image = uiimage(app.UIFigure);
            app.BounLogo_Image.Position = [9 362 100 100];
            app.BounLogo_Image.ImageSource = 'Boun_Logo.png';

            % Create PlayOrigSound_Button
            app.PlayOrigSound_Button = uibutton(app.UIFigure, 'push');
            app.PlayOrigSound_Button.ButtonPushedFcn = createCallbackFcn(app, @PlayOrigSound_ButtonPushed, true);
            app.PlayOrigSound_Button.BackgroundColor = [0 1 1];
            app.PlayOrigSound_Button.FontName = 'Comic Sans MS';
            app.PlayOrigSound_Button.FontSize = 15;
            app.PlayOrigSound_Button.FontWeight = 'bold';
            app.PlayOrigSound_Button.FontColor = [0.149 0.149 0.149];
            app.PlayOrigSound_Button.Position = [65 220 100 71];
            app.PlayOrigSound_Button.Text = {'Play'; 'Original'; 'Speech'};

            % Create PlaySynthSound_Button
            app.PlaySynthSound_Button = uibutton(app.UIFigure, 'push');
            app.PlaySynthSound_Button.ButtonPushedFcn = createCallbackFcn(app, @PlaySynthSound_ButtonPushed, true);
            app.PlaySynthSound_Button.BackgroundColor = [0 1 1];
            app.PlaySynthSound_Button.FontName = 'Comic Sans MS';
            app.PlaySynthSound_Button.FontSize = 15;
            app.PlaySynthSound_Button.FontWeight = 'bold';
            app.PlaySynthSound_Button.FontColor = [0.149 0.149 0.149];
            app.PlaySynthSound_Button.Position = [63 49 103 71];
            app.PlaySynthSound_Button.Text = {'Play'; 'Synthesized'; 'Speech'};

            % Create Save_Button
            app.Save_Button = uibutton(app.UIFigure, 'push');
            app.Save_Button.ButtonPushedFcn = createCallbackFcn(app, @Save_ButtonPushed, true);
            app.Save_Button.BackgroundColor = [0.7176 0.2745 1];
            app.Save_Button.FontName = 'Comic Sans MS';
            app.Save_Button.FontSize = 15;
            app.Save_Button.FontWeight = 'bold';
            app.Save_Button.FontColor = [0.149 0.149 0.149];
            app.Save_Button.Position = [529 174 103 71];
            app.Save_Button.Text = {'Save'; 'Synthesized'; 'Speech'};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LPC_Vocoder

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end