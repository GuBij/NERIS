function class=stab(T114,T8,wspeed)

 GradT = (T114-T8)/106;
 SI = (GradT + 0.0098)/wspeed^2; %*cos(elev*pi/180))^2;

 if (SI ~= 0)
   lambda = log10(abs(SI)*10^6);
 end

 if (wspeed >= 11.5)
   class = 3;
 elseif (SI > 0)
   if (lambda >= 2.75)
	class = 1;
   elseif (lambda > 1.75)
	class = 2;
   else
	class = 3;
   end
 elseif (SI < 0)
   if (lambda <= 2)
	class = 3;
   elseif (lambda < 2.75)
	class = 4;
   elseif (lambda < 3.3)
	class = 5;
   else
	class = 6;
   end
  else
	class = 3;
  end

end
