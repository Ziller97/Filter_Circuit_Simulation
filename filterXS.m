function [p_in, p_out] = filterXS(omega_in, dt, t_max, trans_cutoff)

% Setting up the time domain for analysis. This needs to be done in order
% to compute the input, output signals for desired 'lengths' or data sets.
t = 0 : dt : t_max;

% The input signal has been taken to be sinusoidal curve to represent
% periodic motion that can be useful to study for the FFT analysis and for 
%studying the time domain response and bode plots for the circuit. Basic 
% calculus has been used to calculate dVin.
Vin = sin(2*pi*omega_in*t);
dVin = 2*pi*omega_in*cos(2*pi*omega_in*t);

%This size is important since it will be useful to generate the frequency
%scale and also for performing direct fourier transformation. 
N=size(Vin,2);

%These are the values given for the circuit elements which can be changed
%to cater for different resistors and capacitors. 
R1 = 47000;
R2 = 47000;
R3 = 12000;
C1 = 0.1*10^-6;
C2 = 0.1*10^-6;
C3 = 0.1*10^-6;

%V1,V2, Vout need to be initialised for using for loop to perform matrix
%multiplications for computing the differential function. Here the initial
%conditions have been assumed to be zero for all three nodes, which makes
%sense at time=0, before the input starts, voltages at these nodes should all be zero. 
V1 = zeros(length(t),1);
V2 = zeros(length(t),1);
Vout = zeros(length(t),1);

%Simple Row * Coloumn for matrix calculations have been reduced into a for
%loop to comput all V1, V2 and Vout. The count condition was decided to get
%enough values as the number of values in t, but -1 was added to avoid
%error in exceeding matrix indices when using count+1. 
for count = 1:length(t)-1
    dV1 = (-1/(C2*R3)*V1(count))+(1/(C2*R2)*V2(count))-(1/(C2*R2)*Vout(count))+dVin(count);
    dV2 = (-1/C1)*((1/R1)+(1/R2))*V2(count)+(1/(C1*R2))*Vout(count)+(1/(C1*R1))*(Vin(count));
    dVout = (-1/(C2*R3))*V1(count)+(1/(R2))*((1/C2)+(1/C3))*V2(count)-(1/(R2))*((1/C2)+(1/C3))*Vout(count)+dVin(count);

    V1(count+1) = V1(count) + dt*dV1;
    V2(count+1) = V2(count) + dt*dV2;
    Vout(count+1) = Vout(count) + dt*dVout;
end

% Frequency Scale was established here which can be used for generating
% suitable frequency analysis plots
max_freq = 1000;
df = 1/(dt*N);
freq = 1:df:max_freq;
M = size(freq,2);

%This is where the FFT calculations are done using matlab algorithm. 
%The fft_in is divided into as many values as in Vin. 
fft_in = abs(fft(Vin))/length(Vin);

%The fft out values have been limited from trans_cutoff so that transient
%time response can be skipped and peak value errors can be avoided. 
fft_out = abs(fft(Vout(trans_cutoff:end))/length(Vin));

%Max function is used to find maximum value which will be useful for
%simulating bode plots for circuit analysis. p_in has been restricted in
%length to perform suitable matrix operation later while calculating gain,
%hence must be made to same length.
p_in = max(fft_in(1:length(Vin)));
p_out = max(fft_out);

% plot(t(t<50),Vout(t<50),'Marker','+','MarkerEdgeColor','blue')
% title('V out : Scenario 2')
% xlabel('Time (s)')
% ylabel('V(out) (V)')
end