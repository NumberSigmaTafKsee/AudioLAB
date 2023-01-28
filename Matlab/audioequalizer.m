classdef audio_equalizer_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AudioEqualizerUIFigure        matlab.ui.Figure
        Panel_4                       matlab.ui.container.Panel
        Panel                         matlab.ui.container.Panel
        BrowseButton                  matlab.ui.control.Button
        FilepathEditFieldLabel        matlab.ui.control.Label
        FilepathEditField             matlab.ui.control.EditField
        Panel_3                       matlab.ui.container.Panel
        TypeoffilterButtonGroup       matlab.ui.container.ButtonGroup
        FIRButton                     matlab.ui.control.RadioButton
        IIRButton                     matlab.ui.control.RadioButton
        GOButton                      matlab.ui.control.Button
        FiltertoplotDropDownLabel     matlab.ui.control.Label
        FiltertoplotDropDown          matlab.ui.control.DropDown
        SamplingrateEditField_2Label  matlab.ui.control.Label
        SamplingrateEditField         matlab.ui.control.EditField
        Panel_2                       matlab.ui.container.Panel
        FrequencyBandGainsLabel       matlab.ui.control.Label
        Panel_6                       matlab.ui.container.Panel
        HzSliderLabel                 matlab.ui.control.Label
        gain170                       matlab.ui.control.Slider
        HzSlider_5Label               matlab.ui.control.Label
        gain310                       matlab.ui.control.Slider
        HzLabel                       matlab.ui.control.Label
        gain600                       matlab.ui.control.Slider
        HzLabel_2                     matlab.ui.control.Label
        gain1k                        matlab.ui.control.Slider
        kHzSliderLabel                matlab.ui.control.Label
        gain3k                        matlab.ui.control.Slider
        kHzLabel                      matlab.ui.control.Label
        gain6k                        matlab.ui.control.Slider
        kHzLabel_2                    matlab.ui.control.Label
        gain12k                       matlab.ui.control.Slider
        kHzLabel_3                    matlab.ui.control.Label
        gain14k                       matlab.ui.control.Slider
        kHzLabel_4                    matlab.ui.control.Label
        gain16k                       matlab.ui.control.Slider
        OKButton                      matlab.ui.control.Button
        playButton                    matlab.ui.control.Button
        resumeButton                  matlab.ui.control.Button
        stopButton                    matlab.ui.control.Button
        RESETButton                   matlab.ui.control.Button
        dBLabel                       matlab.ui.control.Label
        dBLabel_2                     matlab.ui.control.Label
        Panel_7                       matlab.ui.container.Panel
        UIAxes2                       matlab.ui.control.UIAxes
        PlotButton                    matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
        Panel_8                       matlab.ui.container.Panel
        UIAxes4                       matlab.ui.control.UIAxes
        PlotSpectrumButton            matlab.ui.control.Button
        UIAxes3                       matlab.ui.control.UIAxes
        FrequencydomainLabel          matlab.ui.control.Label
        TIMEDOMAINLabel               matlab.ui.control.Label
        SaveButton                    matlab.ui.control.Button
        Panel_5                       matlab.ui.container.Panel
        VolumeKnobLabel               matlab.ui.control.Label
        VolumeKnob                    matlab.ui.control.Knob
    end


    properties (Access = private)
        playsound % Description
        y
        fs
        t
        Ns
        fm
        bandgains = db2mag(0)*ones(1,9);
        ycomp = 0;
        filterfs = 48000;
        ycompfs
        newsamplingrate

    end

    methods (Access = private)


        function readFile(app)
            [app.y,app.fs] = audioread(app.FilepathEditField.Value);
            app.Ns = length(app.y);
            app.ycompfs = app.fs;
            app.t = linspace(0, app.Ns/app.fs, app.Ns);
            app.fm = app.fs/2;

        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            movegui(app.AudioEqualizerUIFigure, 'center')
        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            [FileName,FilePath]=uigetfile({'*.wav'});
            fullPath = [FilePath FileName];
            app.FilepathEditField.Value = fullPath;
            readFile(app);
        end

        % Button pushed function: playButton
        function playButtonPushed(app, event)
            global player
                if app.ycomp == 0
                    player = audioplayer(app.y*app.VolumeKnob.Value/100,app.fs);
                    %app.playsound = player;
                else
                    player = audioplayer(app.ycomp*app.VolumeKnob.Value/100,app.ycompfs);
                    %app.playsound = player;
                end
            play(player);
        end

        % Button pushed function: resumeButton
        function resumeButtonPushed(app, event)
            global player
            resume(player);
        end

        % Button pushed function: GOButton
        function GOButtonPushed(app, event)
            freq = [0,170,310,600,1000,3000,6000,12000,14000,16000];
            for x = 1:9
                 title1 = ['Gain and phase response of ', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz filter'];
                 title2 = ['Impulse response of ', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz filter'];
                 title3 = ['Step response of ', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz filter'];
                 title4 = ['poles/zeros of ', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz filter'];
                 title5 = ['Sound signal with (', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz) IIR filter applied in time domain'];
                 title6 = ['Sound signal with (', num2str(freq(x)), ' - ' , num2str(freq(x+1)), ' Hz) FIR filter applied in time domain'];
                 %IIR FILTERS
                 if app.TypeoffilterButtonGroup.SelectedObject == app.IIRButton
                     if app.FiltertoplotDropDown.Value == num2str(x)
                        if app.FiltertoplotDropDown.Value == '1'
                            [num,den] = butter(4,170/(app.filterfs/2));%Filter design in num denom form
                        else
                            [num,den] = butter(4,[freq(x) freq(x+1)]/(app.filterfs/2),'bandpass');
                        end
                        figure;
                        freqz(num,den); title(title1)
                        figure;
                        subplot(2,2,1)
                        impz(num,den); title(title2)
                        subplot(2,2,2)
                        stepz(num,den); title(title3)
                        [z,p] = tf2zpk(num,den);
                        subplot(2,2,[3,4])
                        zplane(z,p); title(title4)

                        filteredSignal = filter(num,den,app.y); %filtered signal in time domain
                        figure;
                        subplot(2,1,1)
                        plot(app.t,filteredSignal)
                        %hold on
                        %plot(app.t,filteredSignal)
                        %legend('Input Data','Filtered Data')
                        title(title5)
                        xlabel('Time (s)'); ylabel('Amplitude');
                        subplot(2,1,2)
                        fmag = abs(fftshift(fft(filteredSignal))/app.fs);
                        L = length(filteredSignal);
                        fvec = linspace(-app.fs/2,app.fs/2,L);
                        plot(fvec,fmag); %filtered signal in frequency domain
                        title('Magnitude of filtered signal in frequency domain'); xlabel('Frequency (Hz)'); ylabel('Magnitude')
                        break
                     end
                %FIR FILTERS
                 else
                     if app.FiltertoplotDropDown.Value == num2str(x)
                        if app.FiltertoplotDropDown.Value == '1'
                            H = fir1(350,170/(app.filterfs/2),'low');%Filter design in num denom form
                        else
                            H = fir1(350,[freq(x) freq(x+1)]/(app.filterfs/2),'bandpass');
                        end
                        figure;
                        freqz(H); title(title1)
                        figure; subplot(2,2,1)
                        impz(H); title(title2)
                        subplot(2,2,2)
                        stepz(H); title(title3)
                        [z,p] = tf2zpk(H,1);
                        subplot(2,2,[3,4])
                        zplane(z,p); title(title4)

                        filteredSignal = filter(H,1,app.y); %filtered signal in time domain
                        figure; subplot(2,1,1)
                        plot(app.t,filteredSignal)
                        title(title6)
                        xlabel('Time (s)'); ylabel('Amplitude');
                        subplot(2,1,2)
                        fmag = abs(fftshift(fft(filteredSignal))/app.fs);
                        L = length(filteredSignal);
                        fvec = linspace(-app.fs/2,app.fs/2,L);
                        plot(fvec,fmag); %filtered signal in frequency domain
                        title('Magnitude of filtered signal in frequency domain'); xlabel('Frequency (Hz)'); ylabel('Magnitude')
                        break
                    end
                end
            end
        end

        % Button pushed function: OKButton
        function OKButtonPushed(app, event)
            gains = [app.gain170.Value app.gain310.Value app.gain600.Value app.gain1k.Value app.gain3k.Value app.gain6k.Value app.gain12k.Value app.gain14k.Value app.gain16k.Value];
            app.bandgains = db2mag(gains);

        end

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
            app.newsamplingrate = app.SamplingrateEditField.Value;
            freq = [0,170,310,600,1000,3000,6000,12000,14000,16000];
            app.ycomp = 0;
            app.filterfs = 48000;
            if app.TypeoffilterButtonGroup.SelectedObject == app.FIRButton
                for x = 1:9
                    if x == 1
                        H = fir1(350,170/(app.filterfs/2));%Filter design in num denom form
                    else
                        H = fir1(350,[freq(x) freq(x+1)]/(app.filterfs/2),'bandpass');
                    end
                    filteredSignal = filter(H,1,app.y); %filtered signal in time domain
                    app.ycomp = app.ycomp + (filteredSignal*app.bandgains(x));
                end
            else
                for x = 1:9
                    if x == 1
                       [num,den] = butter(4,170/(app.filterfs/2));%Filter design in num denom form
                    else
                       [num,den] = butter(4,[freq(x) freq(x+1)]/(app.filterfs/2),'bandpass');
                    end
                    filteredSignal = filter(num,den,app.y); %filtered signal in time domain
                    app.ycomp = app.ycomp + (filteredSignal*app.bandgains(x));
                end
            end
            if strcmpi(app.newsamplingrate,'Default') == 1
                app.ycompfs = app.fs;
            else
                app.ycompfs = str2double(app.newsamplingrate);
            end
            app.ycomp = resample(app.ycomp,app.ycompfs,app.fs);
            plot(app.UIAxes,app.t,app.y)
            tcomp = linspace(0,length(app.ycomp)/app.ycompfs,length(app.ycomp));
            plot(app.UIAxes2,tcomp,app.ycomp)
        end

        % Button pushed function: PlotSpectrumButton
        function PlotSpectrumButtonPushed(app, event)
            app.newsamplingrate = app.SamplingrateEditField.Value;
            %fmag = abs(fftshift(fft(app.y))/app.fs);
            L1 = length(app.y);
            fmag = abs(fftshift(fft(app.y))/L1);
            fvec1 = linspace(-app.fs/2,app.fs/2,L1);
            plot(app.UIAxes4,fvec1,fmag);

            if strcmpi(app.newsamplingrate,'Default') == 1
                app.ycompfs = app.fs;
            else
                app.ycompfs = str2double(app.newsamplingrate);
            end
            app.ycomp = resample(app.ycomp,app.ycompfs,app.fs);
            %fcomp = abs(fftshift(fft(app.ycomp))/app.ycompfs);
            L2 = length(app.ycomp);
            fcomp = abs(fftshift(fft(app.ycomp))/L2);
            fvec2 = linspace(-app.ycompfs/2,app.ycompfs/2,L2);
            plot(app.UIAxes3,fvec2,fcomp);
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            filename = 'D:\Downloads\FILTEREDSOUND.wav';
            audiowrite(filename,app.ycomp,app.ycompfs);
        end

        % Button pushed function: stopButton
        function stopButtonPushed(app, event)
            global player
            pause(player)
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
            app.bandgains = db2mag(0)*ones(1,9);
            app.gain170.Value = 0;
            app.gain310.Value = 0;
            app.gain600.Value = 0;
            app.gain1k.Value = 0;
            app.gain3k.Value = 0;
            app.gain6k.Value = 0;
            app.gain12k.Value = 0;
            app.gain14k.Value = 0;
            app.gain16k.Value = 0;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AudioEqualizerUIFigure and hide until all components are created
            app.AudioEqualizerUIFigure = uifigure('Visible', 'off');
            app.AudioEqualizerUIFigure.Color = [0.0784 0.0784 0.0784];
            app.AudioEqualizerUIFigure.Colormap = [0.2431 0.149 0.6588;0.251 0.1647 0.7059;0.2588 0.1804 0.7529;0.2627 0.1961 0.7961;0.2706 0.2157 0.8353;0.2745 0.2353 0.8706;0.2784 0.2549 0.898;0.2784 0.2784 0.9216;0.2824 0.302 0.9412;0.2824 0.3216 0.9569;0.2784 0.3451 0.9725;0.2745 0.3686 0.9843;0.2706 0.3882 0.9922;0.2588 0.4118 0.9961;0.2431 0.4353 1;0.2196 0.4588 0.9961;0.1961 0.4863 0.9882;0.1843 0.5059 0.9804;0.1804 0.5294 0.9686;0.1765 0.549 0.9529;0.1686 0.5686 0.9373;0.1529 0.5922 0.9216;0.1451 0.6078 0.9098;0.1373 0.6275 0.898;0.1255 0.6471 0.8902;0.1098 0.6627 0.8745;0.0941 0.6784 0.8588;0.0706 0.6941 0.8392;0.0314 0.7098 0.8157;0.0039 0.7216 0.7922;0.0078 0.7294 0.7647;0.0431 0.7412 0.7412;0.098 0.749 0.7137;0.1412 0.7569 0.6824;0.1725 0.7686 0.6549;0.1922 0.7765 0.6235;0.2157 0.7843 0.5922;0.2471 0.7922 0.5569;0.2902 0.7961 0.5176;0.3412 0.8 0.4784;0.3922 0.8039 0.4353;0.4471 0.8039 0.3922;0.5059 0.8 0.349;0.5608 0.7961 0.3059;0.6157 0.7882 0.2627;0.6706 0.7804 0.2235;0.7255 0.7686 0.1922;0.7725 0.7608 0.1647;0.8196 0.749 0.1529;0.8627 0.7412 0.1608;0.902 0.7333 0.1765;0.9412 0.7294 0.2118;0.9725 0.7294 0.2392;0.9961 0.7451 0.2353;0.9961 0.7647 0.2196;0.9961 0.7882 0.2039;0.9882 0.8118 0.1882;0.9804 0.8392 0.1765;0.9686 0.8627 0.1647;0.9608 0.8902 0.1529;0.9608 0.9137 0.1412;0.9647 0.9373 0.1255;0.9686 0.9608 0.1059;0.9765 0.9843 0.0824];
            app.AudioEqualizerUIFigure.Position = [100 100 1134 708];
            app.AudioEqualizerUIFigure.Name = 'Audio Equalizer';
            app.AudioEqualizerUIFigure.Scrollable = 'on';

            % Create Panel_4
            app.Panel_4 = uipanel(app.AudioEqualizerUIFigure);
            app.Panel_4.BorderType = 'none';
            app.Panel_4.BackgroundColor = [0.0784 0.0784 0.0784];
            app.Panel_4.Position = [10 11 1112 686];

            % Create Panel
            app.Panel = uipanel(app.Panel_4);
            app.Panel.BorderType = 'none';
            app.Panel.BackgroundColor = [0.302 0.302 0.3608];
            app.Panel.Position = [20 572 427 80];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.Panel, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.BackgroundColor = [0.2627 0.2941 0.3686];
            app.BrowseButton.FontName = 'Cambria Math';
            app.BrowseButton.FontWeight = 'bold';
            app.BrowseButton.FontColor = [1 1 1];
            app.BrowseButton.Position = [340 27 64 25];
            app.BrowseButton.Text = 'Browse';

            % Create FilepathEditFieldLabel
            app.FilepathEditFieldLabel = uilabel(app.Panel);
            app.FilepathEditFieldLabel.HorizontalAlignment = 'right';
            app.FilepathEditFieldLabel.FontName = 'Calibri';
            app.FilepathEditFieldLabel.FontSize = 14;
            app.FilepathEditFieldLabel.FontWeight = 'bold';
            app.FilepathEditFieldLabel.FontColor = [1 1 1];
            app.FilepathEditFieldLabel.Position = [6 28 62 22];
            app.FilepathEditFieldLabel.Text = 'File path: ';

            % Create FilepathEditField
            app.FilepathEditField = uieditfield(app.Panel, 'text');
            app.FilepathEditField.FontColor = [1 1 1];
            app.FilepathEditField.BackgroundColor = [0.2392 0.2902 0.349];
            app.FilepathEditField.Position = [75 28 226 22];

            % Create Panel_3
            app.Panel_3 = uipanel(app.Panel_4);
            app.Panel_3.BorderType = 'none';
            app.Panel_3.BackgroundColor = [0.302 0.302 0.3608];
            app.Panel_3.Position = [606 572 486 80];

            % Create TypeoffilterButtonGroup
            app.TypeoffilterButtonGroup = uibuttongroup(app.Panel_3);
            app.TypeoffilterButtonGroup.ForegroundColor = [1 1 1];
            app.TypeoffilterButtonGroup.Title = 'Type of filter: ';
            app.TypeoffilterButtonGroup.BackgroundColor = [0.2392 0.2902 0.349];
            app.TypeoffilterButtonGroup.FontName = 'Arial';
            app.TypeoffilterButtonGroup.FontSize = 14;
            app.TypeoffilterButtonGroup.Position = [210 16 176 50];

            % Create FIRButton
            app.FIRButton = uiradiobutton(app.TypeoffilterButtonGroup);
            app.FIRButton.Text = 'FIR';
            app.FIRButton.FontColor = [1 1 1];
            app.FIRButton.Position = [24 2 41 22];
            app.FIRButton.Value = true;

            % Create IIRButton
            app.IIRButton = uiradiobutton(app.TypeoffilterButtonGroup);
            app.IIRButton.Text = 'IIR';
            app.IIRButton.FontColor = [1 1 1];
            app.IIRButton.Position = [114 2 40 22];

            % Create GOButton
            app.GOButton = uibutton(app.Panel_3, 'push');
            app.GOButton.ButtonPushedFcn = createCallbackFcn(app, @GOButtonPushed, true);
            app.GOButton.BackgroundColor = [0.9608 0.949 0.051];
            app.GOButton.FontName = 'Cambria Math';
            app.GOButton.FontSize = 14;
            app.GOButton.FontWeight = 'bold';
            app.GOButton.Position = [421 29 34 24];
            app.GOButton.Text = 'GO';

            % Create FiltertoplotDropDownLabel
            app.FiltertoplotDropDownLabel = uilabel(app.Panel_3);
            app.FiltertoplotDropDownLabel.HorizontalAlignment = 'right';
            app.FiltertoplotDropDownLabel.FontName = 'Calibri';
            app.FiltertoplotDropDownLabel.FontSize = 14;
            app.FiltertoplotDropDownLabel.FontWeight = 'bold';
            app.FiltertoplotDropDownLabel.FontColor = [1 1 1];
            app.FiltertoplotDropDownLabel.Position = [5 16 84 22];
            app.FiltertoplotDropDownLabel.Text = 'Filter to plot: ';

            % Create FiltertoplotDropDown
            app.FiltertoplotDropDown = uidropdown(app.Panel_3);
            app.FiltertoplotDropDown.Items = {'0 - 170 Hz', '170 - 310 Hz', '310 - 600 Hz', '600 - 1000 Hz', '1 - 3 kHz', '3 - 6 kHz', '6 - 12 kHz', '12 - 14 kHz', '14 - 16  kHz'};
            app.FiltertoplotDropDown.ItemsData = {'1', '2', '3', '4', '5', '6', '7', '8', '9'};
            app.FiltertoplotDropDown.FontName = 'Calibri';
            app.FiltertoplotDropDown.FontSize = 14;
            app.FiltertoplotDropDown.FontColor = [1 1 1];
            app.FiltertoplotDropDown.BackgroundColor = [0.2392 0.2902 0.349];
            app.FiltertoplotDropDown.Position = [93 16 103 22];
            app.FiltertoplotDropDown.Value = '1';

            % Create SamplingrateEditField_2Label
            app.SamplingrateEditField_2Label = uilabel(app.Panel_3);
            app.SamplingrateEditField_2Label.HorizontalAlignment = 'right';
            app.SamplingrateEditField_2Label.FontName = 'Calibri';
            app.SamplingrateEditField_2Label.FontSize = 14;
            app.SamplingrateEditField_2Label.FontWeight = 'bold';
            app.SamplingrateEditField_2Label.FontColor = [1 1 1];
            app.SamplingrateEditField_2Label.Position = [5 44 89 22];
            app.SamplingrateEditField_2Label.Text = 'Sampling rate:';

            % Create SamplingrateEditField
            app.SamplingrateEditField = uieditfield(app.Panel_3, 'text');
            app.SamplingrateEditField.FontColor = [1 1 1];
            app.SamplingrateEditField.BackgroundColor = [0.2392 0.2902 0.349];
            app.SamplingrateEditField.Position = [102 44 80 22];
            app.SamplingrateEditField.Value = 'Default';

            % Create Panel_2
            app.Panel_2 = uipanel(app.Panel_4);
            app.Panel_2.ForegroundColor = [0.9412 0.9412 0.9412];
            app.Panel_2.BorderType = 'none';
            app.Panel_2.BackgroundColor = [0.2706 0.2902 0.349];
            app.Panel_2.Position = [20 21 1072 528];

            % Create FrequencyBandGainsLabel
            app.FrequencyBandGainsLabel = uilabel(app.Panel_2);
            app.FrequencyBandGainsLabel.HorizontalAlignment = 'center';
            app.FrequencyBandGainsLabel.FontName = 'ROG Fonts';
            app.FrequencyBandGainsLabel.FontSize = 14;
            app.FrequencyBandGainsLabel.FontColor = [1 1 1];
            app.FrequencyBandGainsLabel.Position = [282 506 227 22];
            app.FrequencyBandGainsLabel.Text = 'Frequency Band Gains';

            % Create Panel_6
            app.Panel_6 = uipanel(app.Panel_2);
            app.Panel_6.BorderType = 'none';
            app.Panel_6.BackgroundColor = [0.3608 0.3804 0.4588];
            app.Panel_6.Position = [7 271 775 221];

            % Create HzSliderLabel
            app.HzSliderLabel = uilabel(app.Panel_6);
            app.HzSliderLabel.HorizontalAlignment = 'center';
            app.HzSliderLabel.FontSize = 10;
            app.HzSliderLabel.FontWeight = 'bold';
            app.HzSliderLabel.FontColor = [1 1 1];
            app.HzSliderLabel.Position = [17 39 52 22];
            app.HzSliderLabel.Text = '0 - 170 Hz';

            % Create gain170
            app.gain170 = uislider(app.Panel_6);
            app.gain170.Limits = [-12 12];
            app.gain170.Orientation = 'vertical';
            app.gain170.FontColor = [1 1 1];
            app.gain170.Position = [29 69 3 135];

            % Create HzSlider_5Label
            app.HzSlider_5Label = uilabel(app.Panel_6);
            app.HzSlider_5Label.HorizontalAlignment = 'center';
            app.HzSlider_5Label.FontSize = 10;
            app.HzSlider_5Label.FontWeight = 'bold';
            app.HzSlider_5Label.FontColor = [1 1 0];
            app.HzSlider_5Label.Position = [92 40 66 22];
            app.HzSlider_5Label.Text = '170 - 310 Hz ';

            % Create gain310
            app.gain310 = uislider(app.Panel_6);
            app.gain310.Limits = [-12 12];
            app.gain310.Orientation = 'vertical';
            app.gain310.FontColor = [1 1 0];
            app.gain310.Position = [108 69 3 135];

            % Create HzLabel
            app.HzLabel = uilabel(app.Panel_6);
            app.HzLabel.HorizontalAlignment = 'center';
            app.HzLabel.FontSize = 10;
            app.HzLabel.FontWeight = 'bold';
            app.HzLabel.FontColor = [0.9294 0.6941 0.1255];
            app.HzLabel.Position = [177 38 62 22];
            app.HzLabel.Text = '310 - 600 Hz';

            % Create gain600
            app.gain600 = uislider(app.Panel_6);
            app.gain600.Limits = [-12 12];
            app.gain600.Orientation = 'vertical';
            app.gain600.FontColor = [0.9294 0.6941 0.1255];
            app.gain600.Position = [194 69 3 135];

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.Panel_6);
            app.HzLabel_2.HorizontalAlignment = 'center';
            app.HzLabel_2.FontSize = 10;
            app.HzLabel_2.FontWeight = 'bold';
            app.HzLabel_2.FontColor = [1 0.1098 0.1098];
            app.HzLabel_2.Position = [267 38 68 22];
            app.HzLabel_2.Text = '600 - 1000 Hz';

            % Create gain1k
            app.gain1k = uislider(app.Panel_6);
            app.gain1k.Limits = [-12 12];
            app.gain1k.Orientation = 'vertical';
            app.gain1k.FontColor = [1 0.1098 0.1098];
            app.gain1k.Position = [281 69 3 135];

            % Create kHzSliderLabel
            app.kHzSliderLabel = uilabel(app.Panel_6);
            app.kHzSliderLabel.HorizontalAlignment = 'center';
            app.kHzSliderLabel.FontSize = 10;
            app.kHzSliderLabel.FontWeight = 'bold';
            app.kHzSliderLabel.FontColor = [1 0.1098 0.451];
            app.kHzSliderLabel.Position = [350 39 46 22];
            app.kHzSliderLabel.Text = '1 - 3 kHz';

            % Create gain3k
            app.gain3k = uislider(app.Panel_6);
            app.gain3k.Limits = [-12 12];
            app.gain3k.Orientation = 'vertical';
            app.gain3k.FontColor = [1 0.1098 0.451];
            app.gain3k.Position = [357 69 3 135];

            % Create kHzLabel
            app.kHzLabel = uilabel(app.Panel_6);
            app.kHzLabel.HorizontalAlignment = 'center';
            app.kHzLabel.FontSize = 10;
            app.kHzLabel.FontWeight = 'bold';
            app.kHzLabel.FontColor = [0.749 0.251 1];
            app.kHzLabel.Position = [420 38 46 22];
            app.kHzLabel.Text = '3 - 6 kHz';

            % Create gain6k
            app.gain6k = uislider(app.Panel_6);
            app.gain6k.Limits = [-12 12];
            app.gain6k.Orientation = 'vertical';
            app.gain6k.FontColor = [0.749 0.251 1];
            app.gain6k.Position = [425 69 3 135];

            % Create kHzLabel_2
            app.kHzLabel_2 = uilabel(app.Panel_6);
            app.kHzLabel_2.HorizontalAlignment = 'center';
            app.kHzLabel_2.FontSize = 10;
            app.kHzLabel_2.FontWeight = 'bold';
            app.kHzLabel_2.FontColor = [0.5098 0.698 1];
            app.kHzLabel_2.Position = [490 39 52 22];
            app.kHzLabel_2.Text = '6 - 12 kHz';

            % Create gain12k
            app.gain12k = uislider(app.Panel_6);
            app.gain12k.Limits = [-12 12];
            app.gain12k.Orientation = 'vertical';
            app.gain12k.FontColor = [0.5098 0.698 1];
            app.gain12k.Position = [502 68 3 135];

            % Create kHzLabel_3
            app.kHzLabel_3 = uilabel(app.Panel_6);
            app.kHzLabel_3.HorizontalAlignment = 'center';
            app.kHzLabel_3.FontSize = 10;
            app.kHzLabel_3.FontWeight = 'bold';
            app.kHzLabel_3.FontColor = [0.4 0.9804 0.9333];
            app.kHzLabel_3.Position = [571 40 57 22];
            app.kHzLabel_3.Text = '12 - 14 kHz';

            % Create gain14k
            app.gain14k = uislider(app.Panel_6);
            app.gain14k.Limits = [-12 12];
            app.gain14k.Orientation = 'vertical';
            app.gain14k.FontColor = [0.4 0.9804 0.9333];
            app.gain14k.Position = [581 68 3 135];

            % Create kHzLabel_4
            app.kHzLabel_4 = uilabel(app.Panel_6);
            app.kHzLabel_4.HorizontalAlignment = 'center';
            app.kHzLabel_4.FontSize = 10;
            app.kHzLabel_4.FontWeight = 'bold';
            app.kHzLabel_4.FontColor = [0.2549 0.9686 0.0392];
            app.kHzLabel_4.Position = [644 38 57 22];
            app.kHzLabel_4.Text = '14 - 16 kHz';

            % Create gain16k
            app.gain16k = uislider(app.Panel_6);
            app.gain16k.Limits = [-12 12];
            app.gain16k.Orientation = 'vertical';
            app.gain16k.FontColor = [0.2549 0.9686 0.0392];
            app.gain16k.Position = [659 69 3 135];

            % Create OKButton
            app.OKButton = uibutton(app.Panel_6, 'push');
            app.OKButton.ButtonPushedFcn = createCallbackFcn(app, @OKButtonPushed, true);
            app.OKButton.BackgroundColor = [1 1 0.0667];
            app.OKButton.FontWeight = 'bold';
            app.OKButton.Position = [342 13 48 26];
            app.OKButton.Text = 'OK';

            % Create playButton
            app.playButton = uibutton(app.Panel_6, 'push');
            app.playButton.ButtonPushedFcn = createCallbackFcn(app, @playButtonPushed, true);
            app.playButton.Icon = 'images.png';
            app.playButton.IconAlignment = 'center';
            app.playButton.BackgroundColor = [0 0 0];
            app.playButton.Position = [726 86 50 50];
            app.playButton.Text = '';

            % Create resumeButton
            app.resumeButton = uibutton(app.Panel_6, 'push');
            app.resumeButton.ButtonPushedFcn = createCallbackFcn(app, @resumeButtonPushed, true);
            app.resumeButton.Icon = 'png-clipart-computer-icons-linkedin-desktop-pause-button-miscellaneous-logo.png';
            app.resumeButton.IconAlignment = 'center';
            app.resumeButton.BackgroundColor = [1 0.9882 0.9882];
            app.resumeButton.Position = [726 151 50 50];
            app.resumeButton.Text = '';

            % Create stopButton
            app.stopButton = uibutton(app.Panel_6, 'push');
            app.stopButton.ButtonPushedFcn = createCallbackFcn(app, @stopButtonPushed, true);
            app.stopButton.Icon = 'stop-button-png-icon-free-download-364870-1.png';
            app.stopButton.IconAlignment = 'center';
            app.stopButton.BackgroundColor = [1 1 1];
            app.stopButton.Position = [726 20 51 50];
            app.stopButton.Text = '';

            % Create RESETButton
            app.RESETButton = uibutton(app.Panel_6, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.BackgroundColor = [1 1 0.0667];
            app.RESETButton.FontWeight = 'bold';
            app.RESETButton.Position = [403 13 55 26];
            app.RESETButton.Text = 'RESET';

            % Create dBLabel
            app.dBLabel = uilabel(app.Panel_6);
            app.dBLabel.FontColor = [1 1 1];
            app.dBLabel.Position = [69 193 25 18];
            app.dBLabel.Text = 'dB';

            % Create dBLabel_2
            app.dBLabel_2 = uilabel(app.Panel_6);
            app.dBLabel_2.FontColor = [1 1 1];
            app.dBLabel_2.Position = [69 48 25 39];
            app.dBLabel_2.Text = ' dB';

            % Create Panel_7
            app.Panel_7 = uipanel(app.Panel_2);
            app.Panel_7.BorderType = 'none';
            app.Panel_7.BackgroundColor = [0.3608 0.3804 0.4588];
            app.Panel_7.Position = [805 46 260 446];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.Panel_7);
            title(app.UIAxes2, 'Composite Signal')
            xlabel(app.UIAxes2, 'Time (s)')
            ylabel(app.UIAxes2, 'Amplitude')
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Position = [6 71 251 156];

            % Create PlotButton
            app.PlotButton = uibutton(app.Panel_7, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.IconAlignment = 'center';
            app.PlotButton.BackgroundColor = [0.9608 0.949 0.051];
            app.PlotButton.FontWeight = 'bold';
            app.PlotButton.Position = [82 22 100 22];
            app.PlotButton.Text = 'Plot';

            % Create UIAxes
            app.UIAxes = uiaxes(app.Panel_7);
            title(app.UIAxes, 'Original Signal')
            xlabel(app.UIAxes, 'Time (s)')
            ylabel(app.UIAxes, 'Amplitude')
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [6 263 251 156];

            % Create Panel_8
            app.Panel_8 = uipanel(app.Panel_2);
            app.Panel_8.BorderType = 'none';
            app.Panel_8.BackgroundColor = [0.3608 0.3804 0.4588];
            app.Panel_8.Position = [7 44 775 221];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.Panel_8);
            title(app.UIAxes4, 'Original  Signal')
            xlabel(app.UIAxes4, 'Frequency (Hz)')
            ylabel(app.UIAxes4, 'Magnitude')
            app.UIAxes4.XGrid = 'on';
            app.UIAxes4.YGrid = 'on';
            app.UIAxes4.Position = [6 24 321 170];

            % Create PlotSpectrumButton
            app.PlotSpectrumButton = uibutton(app.Panel_8, 'push');
            app.PlotSpectrumButton.ButtonPushedFcn = createCallbackFcn(app, @PlotSpectrumButtonPushed, true);
            app.PlotSpectrumButton.BackgroundColor = [0.9608 0.949 0.051];
            app.PlotSpectrumButton.FontWeight = 'bold';
            app.PlotSpectrumButton.Position = [339 3 100 22];
            app.PlotSpectrumButton.Text = 'Plot Spectrum';

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.Panel_8);
            title(app.UIAxes3, 'Composite Signal')
            xlabel(app.UIAxes3, 'Frequency (Hz)')
            ylabel(app.UIAxes3, 'Magnitude')
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.Position = [457 24 305 170];

            % Create FrequencydomainLabel
            app.FrequencydomainLabel = uilabel(app.Panel_8);
            app.FrequencydomainLabel.HorizontalAlignment = 'center';
            app.FrequencydomainLabel.FontName = 'ROG Fonts';
            app.FrequencydomainLabel.FontSize = 14;
            app.FrequencydomainLabel.FontColor = [1 1 1];
            app.FrequencydomainLabel.Position = [297 200 185 22];
            app.FrequencydomainLabel.Text = 'Frequency domain';

            % Create TIMEDOMAINLabel
            app.TIMEDOMAINLabel = uilabel(app.Panel_2);
            app.TIMEDOMAINLabel.HorizontalAlignment = 'center';
            app.TIMEDOMAINLabel.FontName = 'ROG Fonts';
            app.TIMEDOMAINLabel.FontSize = 14;
            app.TIMEDOMAINLabel.FontColor = [1 1 1];
            app.TIMEDOMAINLabel.Position = [877 506 118 22];
            app.TIMEDOMAINLabel.Text = 'TIME DOMAIN';

            % Create SaveButton
            app.SaveButton = uibutton(app.Panel_2, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.BackgroundColor = [0.1098 0.9294 0.1373];
            app.SaveButton.FontSize = 16;
            app.SaveButton.FontWeight = 'bold';
            app.SaveButton.Position = [494 10 100 27];
            app.SaveButton.Text = 'Save';

            % Create Panel_5
            app.Panel_5 = uipanel(app.Panel_4);
            app.Panel_5.BorderType = 'none';
            app.Panel_5.BackgroundColor = [0.0784 0.0784 0.0784];
            app.Panel_5.Position = [470 556 114 107];

            % Create VolumeKnobLabel
            app.VolumeKnobLabel = uilabel(app.Panel_5);
            app.VolumeKnobLabel.HorizontalAlignment = 'center';
            app.VolumeKnobLabel.FontColor = [1 1 1];
            app.VolumeKnobLabel.Position = [34 3 46 22];
            app.VolumeKnobLabel.Text = 'Volume';

            % Create VolumeKnob
            app.VolumeKnob = uiknob(app.Panel_5, 'continuous');
            app.VolumeKnob.FontColor = [0.9294 0.6941 0.1255];
            app.VolumeKnob.Position = [36 34 43 43];
            app.VolumeKnob.Value = 100;

            % Show the figure after all components are created
            app.AudioEqualizerUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = audio_equalizer_GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.AudioEqualizerUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.AudioEqualizerUIFigure)
        end
    end
end

