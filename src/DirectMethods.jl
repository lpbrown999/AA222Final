
#7.5 in book
function nelder_mead(f, S; n=100, α=1.0, β=2.0, γ=0.5)
	#f -> function 
	#S -> Initial simplex 
	#n -> allowed evaluations

	nevals, y_arr = length(S), f.(S)
	while nevals < n
		p = sortperm(y_arr)				#Sorted indexes based on funciton values
		S, y_arr = S[p], y_arr[p]
		xl, yl = S[1], y_arr[1]			#location, lowest value
		xh, yh = S[end], y_arr[end]		#location, highest value
		xs, ys = S[end-1], y_arr[end-1]	#location, second highest
		xm = Statistics.mean(S[1:end-1	])			#Centroid of the simplex
		xr = xm + α*(xm-xh)				#reflection point 
		yr = f(xr);						nevals+=1; 		if (nevals == n) return S, S[argmin(y_arr)] end;

		if yr < yl
			xe = xm + β*(xr-xm)
			ye = f(xe);					nevals+=1; 		if (nevals == n) return S, S[argmin(y_arr)] end;
			S[end], y_arr[end] = ye < yr ? (xe, ye) : (xr, yr)
		elseif yr > ys
			if yr <= yh 
				xh, yh, S[end], y_arr[end] = xr, yr, xr, yr
			end
			xc = xm + γ*(xh - xm)
			yc = f(xc);					nevals+=1; 		if (nevals == n) return S, S[argmin(y_arr)] end;
			if yc > yh
				for i in 2:length(y_arr)
					S[i] = (S[i] + xl)/2
					y_arr[i] = f(S[i]);	nevals+=1; 		if (nevals == n) return S, S[argmin(y_arr)] end;
				end
			else
				S[end], y_arr[end] = xc, yc
			end
		else
			S[end], y_arr[end] = xr, yr
		end
	end
	return S, S[argmin(y_arr)]
end 


#Helpers

#Log indicator 
function log_indicator(x)
	if x <= 0
		return -2^60
	else
		return log(x)
	end
end

# Log barrier for inequality constraints
# fi(x) <= 0. Pass in fi already evaluated
function log_barrier(fi)
	return -(sum(log_indicator.(fi)))
end

#Quadratic penalty for not being hi = 0
function quad_penalty_equality(hi)
	return sum(hi.^2)
end