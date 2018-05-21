figure(2); clf
df_meas = gradient(Fs, deltas);
plot(1:(nn-2), DFs(1:(nn-2)));
hold on
plot(1:(nn-2), df_meas(1:(nn-2)));
legend('Adjoint', 'Calculated');