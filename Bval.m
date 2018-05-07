function [b, b_err, a, R2, X2, Mbins, Ncum, Nbin ]=Bval(M_list, Mc, dM)
  % Bval(M_list, Mc, dM)
  %
  % Function to determine the best fit to the Gutenberg-Richter 
  % frequency-magnitude distribution (GR-FMD).  This code uses the maximum 
  % likelihood method (Aki, 1965) to determine b-value, compensating for 
  % binning error (Marocchi & Sandri, 2003).  Error is determined based on 
  % the method of Shi & Bolt (1982).  Note that the Shi & Bolt method 
  % likely underestimates the true variability in the b-value.
  %
  % Input data:
  % -M_list   -- Magnitudes vector.
  % -Mc       -- Magnitude of completeness.
  % -dM       -- Magnitude bin width and standard error.
  %
  % Output data:
  % -b        -- GR-FMD b-value.
  % -b_err    -- Error in GR-FMD b-value.
  % -a        -- GR-FMD a-value.
  % -R2       -- Goodness-of-fit statistic.
  % -X2       -- Chi-squared statistic.
  % -L        -- Likelihood objective function.
  % -Mbins    -- Magnitude bins, for plotting.
  % -Ncum     -- Cumulative count in GR-FMD, for plotting.
  % -Nbin     -- Frequency count in GR-FMD, for plotting.
  %
  % We kindly ask that if the user finds the use of this code beneficial 
  % that a citation to our research is included in their work:
  %
  % Schultz, R., Atkinson, G., Eaton, D. W., Gu, Y. J., & Kao, H. (2018). 
  % Hydraulic fracturing volume is associated with induced earthquake 
  % productivity in the Duvernay play. Science, 359(6373), 304-308, 
  % doi: 10.1126/science.aao0159.
  %
  % References: 
  %   Aki (1965) Bull. Earthq. Res. Inst. Tokyo Univ. 43,237?239.
  %   Marocchi & Sandri (2003) Ann. Geophys. 46,1271?1282.
  %   Shi & Bolt (1982) Bull. Seismol. Soc. Am. 72(6),1677?1687. 
  %
  %
  % This program is free software: you can redistribute it and/or modify
  % it under the terms of the GNU General Public License as published
  % by the Free Software Foundation, either version 3 of the License, or
  % any later version.
  % This program is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  % GNU General Public License for more details:
  % http://www.gnu.org/licenses/
  
  
  % Ensure that input is a row vector.
  if(iscolumn(M_list))
      M_list=M_list';
  end;
  
  % Check to see there is enough input data points.
  if( length(M_list)<=3 )
      b=NaN; b_err=NaN; R2=NaN; X2=NaN; L=NaN; a=NaN; Mbins=[]; Ncum=[];
      fprintf('Error: Insufficient input data.\n');
      return;
  end;
  
  % Define magnitude bin widths based on the on input error value (dM).
  Mbins=min(M_list)-dM:dM/2:max(M_list)+dM;
  [Nbin,~]=histcounts(M_list, Mbins);
  
  % Get cumulative number of events, but only keep unique bins.
  Ncum=fliplr(cumsum(fliplr(Nbin)));
  diffN=[diff(Ncum),0];
  Ncum=Ncum(diffN~=0);
  Mbins=Mbins(diffN~=0);
  Nbin=Nbin(diffN~=0);
  Ntot=max(Ncum(Mbins>=Mc));
  
  % Define axis for parameter fitting.
  logN=log10(Ncum);
  Y=logN(Mbins>=Mc);
  X=Mbins(Mbins>=Mc);
  
  % Ensure there is still enough data to fit.
  if( length(Y)<=3 )
      b=NaN; b_err=NaN; R2=NaN; X2=NaN; L=NaN; a=NaN; Mbins=[]; Ncum=[];
      fprintf('Error: Insufficient data after magnitude truncation.\n');
      return;
  end;
  
  % Find the mean and standard deviation of the magnitude value.
  temp=[abs(diff(10.^Y)), 10^Y(end)];
  M_avg=sum( X.*temp )/Ntot;
  M_std=sqrt( sum(((X.^2).*temp))/Ntot - M_avg^2 );
  
  % Determine most likely b-value, its error, and the a-value.
  b=1/( log(10)*(M_avg-(Mc-dM/4)) );
  b_err=2.3*M_std*b^2/sqrt(Ntot);
  w=(Nbin(Mbins>=Mc)); w=w/sum(w);
  a=sum((b*X+Y).*w);
  
  % Compute goodness-of-fit statistics (in linear space).
  po=[-b, a];
  Y=10.^Y;
  Yfit=10.^polyval(po,X);
  Ybar=sum(Y.*w);
  SStot=sum(w.*((Y-Ybar).^2));
  SSres=sum(w.*((Y-Yfit).^2));
  R2=1-(SSres/SStot);
  X2=SSres/(b_err);
  
return;
