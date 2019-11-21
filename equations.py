import numpy

def logistic_equation(t, k, N, r) -> float:
	A = (k - N) / N
	E = numpy.exp(-r * t)

	result = k / (1 + (A * E))
	return result


def logistic_equation_integral(t, k, N, r) -> float:
	A = (k - N) / N
	numerator = k * numpy.log(A + numpy.exp(r * t))

	return numerator / r

if __name__ == "__main__":
	pass