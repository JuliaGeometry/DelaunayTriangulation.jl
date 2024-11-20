@inline check_precision(x) = abs(x) ≤ ε(x)

@inline check_absolute_precision(x, y) = check_precision(x - y)

@inline function check_relative_precision(x, y)
    x, y = abs(x), abs(y)
    if x < y
        x, y = y, x
    end
    return !iszero(x) && check_precision(abs(x - y) / x)
end

const RATIO_LB = 0.99
const RATIO_UB = 1.01
@inline check_ratio_precision(x, y) = !iszero(y) && RATIO_LB < abs(x / y) < RATIO_UB