import numpy
from loguru import logger

very_small_number = 1E-6


def logistic_equation(t, k, N, r) -> float:
	A = (k - N) / N
	E = numpy.exp(-r * t)

	result = k / (1 + (A * E))
	return result


def logistic_equation_integral(t, k, N, r) -> float:
	A = (k - N) / N
	numerator = k * numpy.log(A + numpy.exp(r * t))
	result = numerator / r
	if result < 0:
		logger.warning(f"The calculated AUC is negative: {result}. Using a value of {very_small_number} instead.")
		result = very_small_number
	return result


if __name__ == "__main__":
	pass
