function [data,fsamp] = commonPreProc(data)
        %% re-reference data
        fsamp = 1000;
        ampSF = 0.0298; % amp units -> 0.0298 mV
        data = data.*ampSF; % convert to real units
        data=car(data); % common average re-referencing (K. Miller Function)
        data = bandpass(data,[4 98],fsamp); % thresh is -2.5
