// --+ options: json=compute +--
var ffr, unrate, cpi;
varexo e_ffr, e_unrate, e_cpi;

model;

[eqnum='ffr']
     ffr = adl(ffr, 'p_ffr_ffr', [1:3]) + adl(unrate, 'p_ffr_unrate', 1) + adl(cpi, 'p_ffr_cpi', [4]);

[eqnum='unrate']
     unrate = adl(unrate, 'p_ffr_unrate', [4 2 5]) + adl(cpi, 'p_unrate_cpi', 6);

[eqnum='cpi']
     cpi = adl(ffr, 'p_cpi_ffr', 2) + adl(cpi, 'p_cpi_cpi', [2]);

end;

// Must be calibrated after the model block
p_ffr_ffr_lag_1 = 1;
p_ffr_ffr_lag_2 = p_ffr_ffr_lag_1*.5;
p_ffr_ffr_lag_3 = p_ffr_ffr_lag_2*.5;
p_ffr_unrate_lag_1 = 2;
p_ffr_cpi_lag_4 = 3;

// Actual paths for the variables.
ds1 = dseries(randn(30, 3), 1, {'ffr', 'unrate', 'cpi'});

// Baseline paths for the variables.
ds0 = dseries(zeros(30, 3), 1, {'ffr', 'unrate', 'cpi'});

olseqs(ds1, 'eqnum', {'ffr', 'cpi'});

sur(ds1);

surgibbs(ds1, randn(17,17), 1000);

plot_contributions('eqnum', 'ffr', ds1, ds0);