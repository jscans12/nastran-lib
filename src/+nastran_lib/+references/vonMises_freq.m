function vonMises = vonMises_freq(sigmaX,sigmaY,sigmaZ,tauXY,tauYZ,tauZX,PSDF)
%VONMISES_FREQ Calculate frequency-domain Von Mises stress
%
%http://www.vibrationdata.com/tutorials_alt/Sandia_vm_stress.pdf
    
    % Defaults
    if isempty(sigmaZ)
        sigmaZ = zeros(size(sigmaX));
    end
    if isempty(tauYZ)
        tauYZ = zeros(size(sigmaX));
    end
    if isempty(tauZX)
        tauZX = zeros(size(sigmaX));
    end
    if nargin < 7
        PSDF = [];
    end
    
    % Build the A-matrix per Sandia paper
    A = [ 1.0 -0.5 -0.5  0.0  0.0  0.0
         -0.5  1.0 -0.5  0.0  0.0  0.0
         -0.5 -0.5  1.0  0.0  0.0  0.0
          0.0  0.0  0.0  3.0  0.0  0.0
          0.0  0.0  0.0  0.0  3.0  0.0
          0.0  0.0  0.0  0.0  0.0  3.0];
    
    % Build the stress matrix
    sigma = [sigmaX,sigmaY,sigmaZ,tauXY,tauYZ,tauZX]';
    
    % Calculate Von Mises stress
    for i = 1:size(sigma,2)
        vonMises(i) = real(sigma(:,i)'*A*sigma(:,i));
    end
    
    % Multiply by forcing function
    if ~isempty(PSDF)
        vonMises = vonMises .* PSDF';
    end
    
end

