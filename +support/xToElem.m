function elem = xToElem(x, xNodes)

[n, bin] = histc(x, xNodes);

elem = bin;