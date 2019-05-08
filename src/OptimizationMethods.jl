import Statistics
#First order optimization method
function nelder_mead(f, S, n; α=1.0, β=2.0, γ=0.5)
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


## Penalty methods
# c a value supposed to be <= 0.

function quad_penalty(c)
    return sum( max.(0,c).^2 )
end

function count_penalty(c)
    return sum( c .>= 0 )
end

function combined_penalty(c,p1,p2)
    return p1*quad_penalty(c) + p2*count_penalty(c)
end

## Generate a simplex using uniform distribution from a to b
function generate_simplex(ndim::Int64,a::Float64,b::Float64)
    S = [ a.+(b-a)*rand(ndim) for i in 1:ndim+1 ]
    return S
end
#If each dim needs dif range
function generate_simplex(ndim::Int64,a::Vector,b::Vector)
	S = []
	for i in 1:ndim+1
		push!(S,zeros(ndim))
		for j in 1:ndim
			S[i][j] = a[j].+(b[j]-a[j])*rand()
		end
	end
	return S
end