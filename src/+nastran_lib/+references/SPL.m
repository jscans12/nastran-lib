function pointspl = SPL(p_Pa,p_ref)
%SPL Calculate sound pressure level (SPL)

    p_rms = sqrt(mean(p_Pa.^2));
    pointspl = 20*log10(p_rms/p_ref);

end

