function plotRawTimeSeries_burst(R,data,mod)
tvec_obs = R.IntP.tvec_obs;
tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
R.IntP.tvec_obs = tvec_obs;
plot(repmat(R.IntP.tvec_obs,size(data,1),1)',data');
xlabel('Time (s)'); ylabel('Amplitude')

