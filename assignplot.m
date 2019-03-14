loglog(1:0.1:1000,NaN(1,9991));
hold on;

dt=0.0001;
t_max=1000;
trans_cutoff=10;

for frequency = [1:1:500 500:100:1000]
    [p_in, p_out] = filterXS(frequency,dt,t_max,trans_cutoff);
    gain = p_out/p_in;
    plot(frequency,gain,'Marker','o');
    title('Bode Plot : Scenario 2')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude/ Gain (dB)')
    drawnow;
end