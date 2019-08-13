from pypka import Titration
from cli import checkParsedInput

# Assignment of global variables
parameters = checkParsedInput()

Titration(parameters, sites=None)

print('CLI exited successfully')
