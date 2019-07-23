% EOF analysis for annual data

function [COEFF,SCORE,rotscores,latent,tsquare,var_explained] = EOF_calc(EOF_data)

[COEFF,SCORE,latent,tsquare] = princomp(EOF_data);
var_explained = cumsum(latent)./sum(latent)*100;

% COEFF'*COEFF %check for orthogonality

rotatedcoeffs = rotatefactors(COEFF(:,1:6),'Method','varimax');
rotscores = zscore(EOF_data)*rotatedcoeffs;

end



