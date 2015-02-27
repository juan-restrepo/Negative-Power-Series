class PowerSeries:

    """Infinite polynomials, with a bounded precision. It supports
    addition, multiplication, subtraction, division (when the leading
    term is 1 or -1; this can be changed easily, but this way the
    coefficients will always be integers) and multiplication by 
    scalars.
    """

	def __init__(self, coeffs, negcoeffs = [], prec = 10):
		
		M = len(coeffs)
		if M < prec:
			zeros = [0 for i in range(prec-M)]
			self.coeffs = coeffs + zeros
		else:
			self.coeffs = coeffs[:prec]
		
		if all(term == 0 for term in negcoeffs):
			self.negcoeffs = []
		else:
			while negcoeffs[-1] == 0:
				del negcoeffs[-1]
			self.negcoeffs = negcoeffs
		
		self.prec = prec
		

	def __add__(self, other):
		m = min([self.prec, other.prec])
		powersum = []
		for i in range(m):
			powersum.append(self.coeffs[i] + other.coeffs[i])
		lenselfneg = len(self.negcoeffs)
		lenotherneg = len(other.negcoeffs)
		powersumneg = []

		if lenselfneg >= lenotherneg:
			rest = self.negcoeffs[lenotherneg:lenselfneg]
		else:
			rest = other.negcoeffs[lenselfneg:lenotherneg]

		for i in range(min(lenotherneg,lenselfneg)):
			powersumneg.append(self.negcoeffs[i] + other.negcoeffs[i])
		powersumneg += rest

		return PowerSeries(powersum,powersumneg,m)


	def __sub__(self, other):
		m = min([self.prec, other.prec])
		powersub = []
		for i in range(m):
			powersub.append(self.coeffs[i] - other.coeffs[i])
		lenselfneg = len(self.negcoeffs)
		lenotherneg = len(other.negcoeffs)
		powersubneg = []

		if lenselfneg >= lenotherneg:
			rest = self.negcoeffs[lenotherneg:lenselfneg]
		else:
			rest = []
			for i in range(lenselfneg,lenotherneg):
				rest.append(-other.negcoeffs[i])

		for i in range(min(lenotherneg,lenselfneg)):
			powersubneg.append(self.negcoeffs[i] - other.negcoeffs[i])
		powersubneg += rest

		return PowerSeries(powersub,powersubneg,m)

	def __mul__(self, other):
		
		poleself = len(self.negcoeffs)
		poleother = len(other.negcoeffs)
		zeroself = 0
		zeroother = 0
		if poleself == 0:
			newself = self.coeffs + []
			while newself[0] == 0:
				del newself[0]
				zeroself += 1
		else:
			newself = self.negcoeffs[::-1] + self.coeffs	
		if poleother == 0:
			newother = other.coeffs + []
			while newother[0] == 0:
				del newother[0]
				zeroother += 1
		else:
			newother = other.negcoeffs[::-1] + other.coeffs
		
		m = min(len(newself), len(newother))
		powermul = []
		
		for i in range(m):
			coef = 0
			for j in range(i+1):
				coef += newself[j] * newother[i-j]
			powermul.append(coef)

		exponent = zeroself + zeroother - poleself - poleother

		if exponent == 0:
			return PowerSeries(powermul, [] , m)
		elif exponent > 0:
			return PowerSeries([0 for i in range(exponent)] + powermul, [] , m + exponent)
		return PowerSeries(powermul[-exponent:], powermul[-exponent-1::-1], m + exponent)

	def inverse(self):

		'''	Finds the inverse of a power series.
		Only works if first term is 1 or -1.
		'''

		pole = len(self.negcoeffs)
		zero = 0
		coeffs = self.negcoeffs[::-1] + self.coeffs
		while coeffs[0] == 0:
			del coeffs[0]
			zero += 1
			
		assert abs(coeffs[0]) == 1, 

		if pole == 0:
			M = self.prec - zero
		else:
			M = self.prec + pole

		inv = [1/coeffs[0]]
		for i in range(1,M):
			coef = 0
			for j in range(1,i+1):
				coef -= coeffs[j] * inv[i-j]
			coef /= coeffs[0]
			inv.append(coef)

		if pole == 0:
			if zero > 0:
				return PowerSeries(inv[zero:],inv[zero - 1::-1], M - zero)
			return PowerSeries(inv[zero:],[], M)
		zeros = [0 for i in range(pole)]
		inv = zeros + inv
		return PowerSeries(inv, [], M + pole)

	def __div__(self, other):

		return self*other.inverse()

	def exp(self, exponent):

		'''Fast exponentiation algorithm.'''

		if exponent == 0:
			identity = [1]
			for j in range(1,self.prec):
				identity.append(0)
			return PowerSeries(identity, [], self.prec)
		if exponent == 1:
			return self
		if exponent%2 == 0:
			return self.exp(exponent/2)*self.exp(exponent/2)
		return self.exp(exponent/2)*self.exp(exponent/2)*self

	def __eq__(self, other):

		'''Checks whether or not two power series are equal up to the largest possible precision given by them.'''

		prec = min(self.prec, other.prec)
		if prec <= 0:
			return self.negcoeffs[-prec:] == other.negcoeffs[-prec:]
		return self.negcoeffs == other.negcoeffs and self.coeffs[:prec] == other.coeffs[:prec]

	def scalarmult(self, c):

		'''Multiplies the power series by the scalar c.'''

		newcoeffs = [c*coeff for coeff in self.coeffs]
		newnegcoeffs = [c*coeff for coeff in self.negcoeffs]

		return PowerSeries(newcoeffs, newnegcoeffs, self.prec)

	def printpole(self):

		"""Prints only the part with poles
(coefficients with negative exponents).
		"""

		pole = len(self.negcoeffs)
		if pole == 0:
			return "O(1)"
		if self.negcoeffs[-1] == 1:
			string = "z^-%d" % pole
		elif self.negcoeffs[-1] == -1:
			string = "-z^-%d" % pole
		else:
			string = "%dz^-%d" % (self.negcoeffs[-1], pole)

		m = min(0, self.prec)
		for index in range(pole-2,-m-1,-1):
			if self.negcoeffs[index] == 0:
				continue
			elif self.negcoeffs[index] == 1:
				string += " + z^-%d" % (index+1)
			elif self.negcoeffs[index] == -1:
				string += " - z^-%d" % (index+1)
			elif self.negcoeffs[index] > 0:
				string += " + %dz^-%d" % (self.negcoeffs[index], index+1)
			else:
				string += " - %dz^-%d" % (-self.negcoeffs[index], index+1)
		if m < 0:
			return string + " + O(z^%d)" % self.prec
		
		return string + " + O(1)"

	def __str__(self):

		pole = len(self.negcoeffs)
		if pole == 0:
			haspole = False
		else:
			haspole = True
		if haspole:
			if self.negcoeffs[-1] == 1:
				string = "z^-%d" % pole
			elif self.negcoeffs[-1] == -1:
				string = "-z^-%d" % pole
			else:
				string = "%dz^-%d" % (self.negcoeffs[-1], pole)

			m = min(0, self.prec)
			for index in range(pole-2,-m-1,-1):
				if self.negcoeffs[index] == 0:
					continue
				elif self.negcoeffs[index] == 1:
					string += " + z^-%d" % (index+1)
				elif self.negcoeffs[index] == -1:
					string += " - z^-%d" % (index+1)
				elif self.negcoeffs[index] > 0:
					string += " + %dz^-%d" % (self.negcoeffs[index], index+1)
				else:
					string += " - %dz^-%d" % (-self.negcoeffs[index], index+1)
			if m == self.prec:
				if m == 0:
					return string + " + O(1)"
				else:
					return string + " + O(z^%d)" % self.prec
		index = 0
		while(self.coeffs[index] == 0):
			index += 1
			if index == self.prec:
				if haspole:
					return string + " + O(z^%d)" % self.prec
				else:
					return "O(z^%d)" % self.prec
		if index == 0:
			if haspole:
				if self.coeffs[index] > 0:
					string += " + %d" % self.coeffs[index]
				else:
					string += " - %d" % -self.coeffs[index]
			else:
				string = str(self.coeffs[index])
		else:
			if self.coeffs[index] == 1:
				if index == 1:
					if haspole:
						string += " + z"
					else:
						string = "z"
				else:
					if haspole:
						string += " + z^%d" % index
					else:
						string = "z^%d" % index
			elif self.coeffs[index] == -1:
				if index == 1:
					if haspole:
						string += " - z"
					else:
						string = "-z"
				else:
					if haspole:
						string += "- z^%d" % index
					else:
						string = "-z^%d" % index
			else:
				if index == 1:
					if haspole:
						if self.coeffs[1] > 0:
							string += " + %dz" % self.coeffs[1]
						else:
							string += " - %dz" % -self.coeffs[1]
					else:
						string = str(self.coeffs[index])+"z"
				else:
					if haspole:
						if self.coeffs[index] > 0:
							string += " + %dz^%d" % (self.coeffs[index], index)
						else:
							string += " - %dz^%d" % (-self.coeffs[index], index)
					else:
						string = str(self.coeffs[index])+"z^%d" % index
		for i in range(index+1,self.prec):
			if i == 1:
				if self.coeffs[i] > 0:
					if self.coeffs[i] == 1:
						string += " + z"
					else:
						string += " + %dz" % self.coeffs[i]
				elif self.coeffs[i] < 0:
					if self.coeffs[i] == -1:
						string += " - z"
					else:
						string += " - %dz" % -self.coeffs[i]
			else:
				if self.coeffs[i] > 0:
					if self.coeffs[i] == 1:
						string += " + z^%d" % i
					else:
						string += " + %dz^%d" % (self.coeffs[i], i)
				elif self.coeffs[i] < 0:
					if self.coeffs[i] == -1:
						string += " - z^%d" % i
					else:
						string += " - %dz^%d" % (-self.coeffs[i], i)
		string += " + O(z^%d)" % self.prec
		return string


	''' This multiplication method works too, but somehow has less precision than the one chosen above.

	def multiplication(self, other):
		poleself = len(self.negcoeffs)
		poleother = len(other.negcoeffs)
		m = min([self.prec+poleself, other.prec+poleother])
		powermul = []
		newself = self.negcoeffs[::-1] + self.coeffs
		newother = other.negcoeffs[::-1] + other.coeffs
		for i in range(m):
			coef = 0
			for j in range(i+1):
				coef += newself[j] * newother[i-j]
			powermul.append(coef)

		if poleself + poleother > 0:
			return PowerSeries(powermul[poleself + poleother :], powermul[poleself + poleother - 1::-1] ,m-poleself-poleother)
		return PowerSeries(powermul, [] ,m)
		'''

def main():
	vec = [5,2,-2,5,6,2,-2,5,3,6,2,6,3,-6,-7,-21,43]
	negvec = [4,-2,1]
	othervec = [-3,6,5,2,5,7,43,2,5]
	othernegvec = [6,2,3,-6,2]
	vec0 = [0,6,2,6,8,3,2,7,5]
	vec000 = [0,0,0,8,3,6,2,7,3,6,67]
	vec00000 = [0,0,0,0,0,1,2,3,4,5,6,7,8,9]

	f = PowerSeries(vec, negvec, 14)
	# 1/x^3-2/x^2+4/x+5+ 2x -2x^2+ 5x^3+ 6x^4+2x^5 -2x^6+5x^7+3x^8 +6x^9+2x^10+6x^11+3x^12-6x^13-7x^14-21x^15+43x^16
	# 1/z^3-2/z^2+4/z+5+2z-2z^2+5z^3+6z^4+2z^5-2z^6+5z^7+3z^8+6z^9+2z^10+6z^11+3z^12-6z^13+O(z^14)
	g = PowerSeries(othervec, othernegvec, 12)
	# 2/x^5-6/x^4+3/x^3+ 2/x^2+6/x-3+ 6x+ 5x^2+ 2x^3+ 5x^4+7x^5+43x^6 +2x^7+5x^8
	# 2/z^5-6/z^4+3/z^3+2/z^2+6/z-3+6z+5z^2+2z^3+5z^4+7z^5+43z^6+2z^7+5z^8+O(z^9)


	# print "f =", f
	# print "g =", g

	# h = f + g
	# print "f + g =", h.negcoeffs[::-1], h.coeffs, h.prec

	# h = g + f
	# print "g + f =", h.negcoeffs[::-1], h.coeffs, h.prec

	# h = f - g
	# print "f - g =", h.negcoeffs[::-1], h.coeffs, h.prec

	# h = g - f
	# print "g - f =", h.negcoeffs[::-1], h.coeffs, h.prec

	# h = f * g
	# print "f * g =", h

	# h = g * f
	# print "g * f =", h

	zero = PowerSeries(vec0 , [], len(vec0))
	zero3 = PowerSeries(vec000 , [], len(vec000))
	zero5 = PowerSeries(vec00000 , [], len(vec00000))

	# vec0 = 6x+2x^2+6x^3+8x^4+3x^5+2x^6+7x^7+5x^8
	# vec000 = 8x^3+3x^4+6x^5+2x^6+7x^7+3x^8+6x^9+67x^10
	# vec00000 = x^5+2x^6+3x^7+4x^8+5x^9+6x^10+7x^11+8x^12+9x^13

	# vec0 = 6z+2z^2+6z^3+8z^4+3z^5+2z^6+7z^7+5z^8+O(z^9)
	# vec000 = 8z^3+3z^4+6z^5+2z^6+7z^7+3z^8+6z^9+67z^10+O(z^11)
	# vec00000 = z^5+2z^6+3z^7+4z^8+5z^9+6z^10+7z^11+8z^12+9z^13+O(z^14)

	gzero = PowerSeries(vec0 , othernegvec, len(vec0))
	gzero3 = PowerSeries(vec000 , othernegvec, len(vec000))
	gzero5 = PowerSeries(vec00000 , othernegvec, len(vec00000))

	# 2/z^5-6/z^4+3/z^3+2/z^2+6/z+6z+2z^2+6z^3+8z^4+3z^5+2z^6+7z^7+5z^8+O(z^9)
	# 2/z^5-6/z^4+3/z^3+2/z^2+6/z+8z^3+3z^4+6z^5+2z^6+7z^7+3z^8+6z^9+67z^10+O(z^11)
	# 2/z^5-6/z^4+3/z^3+2/z^2+6/z+z^5+2z^6+3z^7+4z^8+5z^9+6z^10+7z^11+8z^12+9z^13+O(z^14)

	# zeros = (zero, zero3, zero5, gzero, gzero3, gzero5)

	# for series in zeros:

	# 	print "series = ", series
	# 	h = f * series
	# 	print "f * series =", h

	# 	h = g * series
	# 	print "g * series =", h

	print "f = ", zero5
	print "1/f = ", zero5.inverse()
	print "f * 1/f = ", zero5*zero5.inverse()
	print "f/f =", zero5/zero5
	# print "f prod 1/f = ", zero5.multiplication(zero5.inverse())
	print "g = ", g
	print "f * g = ", zero5 * g
	# print "f prod g = ", zero5.multiplication(g)

	gg = PowerSeries(othervec, othernegvec, 4)
	print "gg = ", gg
	print "f * gg = ", zero5 * gg
	# print "f prod gg = ", zero5.multiplication(gg)

	print "\n\n"

	print "This is 1/(1-x-x^2), which has the n-th Fibonacci number as the coefficient of z^n."

	print PowerSeries([1,-1,-1], prec = 20).inverse()





if __name__=='__main__':
	main()
