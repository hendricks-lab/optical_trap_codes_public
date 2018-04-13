# optical_trap_codes_public
Processing and analysis of optical trapping data


Step1_TransferFunctionCalculation.m - Calculates the transfer function of the quadrant photodiode signal (V) to sinusoidal excitations of either the optical trap position or the stage position. This code should be used for a sweep of single-harmonic excitations.

Fitting_adam_method_r2.m - This code fits a simple viscoelastic model to the measured transfer function (from the code above) to find the trap stiffness and estimate the viscoelastic properties of the media (see Hendricks et al., 2012, PNAS).

