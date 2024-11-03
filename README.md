The DUSTY model has two parts. Run "cal_0_Ma_dust.ncl" first to tune the model to pre-industrial and derive the scaling constants, namely C_1, C_2, C_3, C_4, and C_5 (see details in the paper). Then run "cal_paleo_dust.ncl," which adapts those scaling constants and applies to all the paleo experiments.

The factorisation analysis is also embedded in the model code.


------
Unzip the dust_emission_fields.zip  to get dust emission fields for all 109 paleo time slices.

The correspondence between the exp name and the simulated paleo time is in name&age.txt
