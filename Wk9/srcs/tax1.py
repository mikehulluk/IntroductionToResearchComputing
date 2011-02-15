



def tax( i ):
  return i * 0.05

def tax( i ):
  if i < 10000:
    return i * 0.05
  else:
    return (i - 10000) * 0.10 + 10000 * 0.05

def tax( i, a ):
  if i < 10000:
    if a > 65:
      return i * 0.04
    else:
      return i * 0.05
  else:
    if a > 65:
      return (i - 10000) * 0.10 + 10000 * 0.04
    else:
      return (i - 10000) * 0.10 + 10000 * 0.05


def tax( i, a, ss ):
  if ss:
    if i < 10000:
      if a > 65:
        return i * 0.03
      else:
        return i * 0.04
    else:
      if a > 65:
        return (i - 10000) * 0.09 + 10000 * 0.03
      else:
        return (i - 10000) * 0.09 + 10000 * 0.04

  else:
    if i < 10000:
      if a > 65:
        return i * 0.04
      else:
        return i * 0.05
    else:
      if a > 65:
        return (i - 10000) * 0.10 + 10000 * 0.04
      else:
        return (i - 10000) * 0.10 + 10000 * 0.05







def getHigherBandTaxRate( age, specialstatus ):
  isRetired = age > 65

  if isRetired and specialstatus: 
    return 0.09

  if isRetired and not specialstatus:
    return 0.10

  if not isRetired and specialstatus:
    return 0.09

  if not isRetired and not specialstatus:
    return 0.10
  

def getLowerBandTaxRate( age, specialstatus):
  isRetired = age > 65

  if isRetired and specialstatus: 
    return 0.03

  if isRetired and not specialstatus:
    return 0.04

  if not isRetired and specialstatus:
    return 0.04

  if not isRetired and not specialstatus:
    return 0.05


def calcTaxToPay(income, age, specialstatus):
  taxBandThreshold = 10000

  if income > taxBandThreshold:
    return taxBandThreshold * getLowerBandTaxRate(age, income) + (income -taxBandThreshold) * getHigherBandTaxRate
  else
    return taxBandThreshold * getLowerBandTaxRate(age, income) 


