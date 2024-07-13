%this function computes the sample entropy for the given data vector 
%with r tolerance and delay

%input: data: moving average of the data points
%          r: tolerance for similarity 
%      delay: shift(1:scale), more shift for higher scale
%output: entropy


function entropy = SampEn(data,r,delay);
n = length(data);
nm   = 0;   %number of matched samples of dimension m
nmp1 = 0;   %number of matched samples of dimension m+1

m = 2;
for i = 1:n-(m+1)*delay;
  for j = i+delay:1:n-m*delay;
    if abs(data(i)-data(j)) < r && abs(data(i+(m-1)*delay)-data(j+(m-1)*delay)) < r; 
      nm = nm + 1;
      if abs(data(i+m*delay)-data(j+m*delay)) < r;
        nmp1 = nmp1+1;
      end
    end
  end
end

%------------------------------------------------------------------------
%nmp1 and nm can be scaled by (n-m+1) and (n-m), but that does not matter 
%as we scale nmp1 and nm both by the same factor.
%------------------------------------------------------------------------
entropy = -log(nmp1/nm);


end
