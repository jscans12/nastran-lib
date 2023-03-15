function srs_aa = srs(data,dt,Q,freq)
%SRS Acceleration SRS using Smallwood ramp invariant digital recursive
%filtering relationship
%
%The highest natural frequency should be less than the Nyquist frequency.
%Ideally, the highest natural frequency should less that one-tenth of
%the sample rate.
%
%srs = [ natural frequency(Hz), peak pos(G), absolute value of peak neg(G)]
%
%INPUTS:
%------
%
%data - input acceleration array - single column - accel (G)
%dt   - constant time step
%Q    - damping quality factor
%freq - array of natural frequencies Hz (pick your own)
%
%--------------------------------------------------------------------------
%
% ver 1.2  May 15, 2014
%
%By Tom Irvine
%Adapted for this toolbox J.Scanlon 
    
    % num_fn = number of natural frequencies
    num_fn = max(size(freq));
    damp = 1/(2*Q);
    
    % The first loop calculates the Smallwood coefficients
    a1 = zeros(num_fn,1);
    a2 = zeros(num_fn,1);
    b1 = zeros(num_fn,1);
    b2 = zeros(num_fn,1);   
    b3 = zeros(num_fn,1);   
    for j = 1:num_fn   
        
        omega = (2*pi*freq(j));

        omegad = (omega*sqrt(1.-damp^2));

        E = (exp(-damp*omega*dt));
        K = (omegad*dt);
        C = (E*cos(K));
        S = (E*sin(K));

        Sp = S/K;

        a1(j) = (2*C);
        a2(j) = (-(E^2));

        b1(j) = (1.-Sp);
        b2(j) = (2.*(Sp-C));
        b3(j) = ((E^2)-Sp);
            
    end
    
    % The second loop calculates the SRS 
    % srs_aa = peak response abs acceleration
    srs_aa = zeros(num_fn,1);
    for j = 1:num_fn

        forward = [ b1(j),  b2(j),  b3(j) ];    
        back    = [     1, -a1(j), -a2(j) ];    

        resp = filter(forward,back,data);

        srs_aa(j) = max(abs(resp));
            
    end
    
end

