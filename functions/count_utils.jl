# Utility functions for conveient math/counting operations

# Repeat elements of a vector x, where each element is repeated the number
# of times specified in the corresponding element of the vector lenghts.
function rep(x, lengths)
    res = similar(x, sum(lengths))
    i = 1
    for xi in 1:length(x)
        val = x[xi]
        for j in 1:lengths[xi]
            res[i] = val
            i += 1
        end
    end
    return res
end

# Find the two integers closest together whos product is num.
function closest_divisors(num)
	b = []
	d = []
	for i in 1:num
		if rem(num,i) == 0
			append!(b,i);
		end
	end
	if rem(size(b,1),2) != 0
		return repeat([median(b)],2)
	else return [b[Int(size(b,1)/2)],b[Int(size(b,1)/2 + 1)]]
	end
end

# Find the two integers who have, at most, a difference of 1 and whos product
# is approximately (if not exactly) num.
function split_int(num)
	[floor(sqrt(num)), ceil(sqrt(num))]
end


# Return boolean indicating whether x is between the min and max of a
# vector range.
function is_between(x,range)
	x >= minimum(range) && x <= maximum(range)
end

# Take square root safely (if negative returns 0)
function sqrt_zero(x)
        sqrt(max(x,0))
end
