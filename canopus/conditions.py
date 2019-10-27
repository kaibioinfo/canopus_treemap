import re
class Condition(object):

	def __init__(self, name, regexp=None, color=None):
		self.name = name
		self.regexpString = regexp if regexp else ".*"+regexp+".*"
		self.regexp = re.compile(regexpString)
		self.color = color
		self.fallbackColor = None

	def coloring(self):
		return self.color if self.color else self.fallbackColor

