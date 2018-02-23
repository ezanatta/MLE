from rcolf import plots_from_prob as pp
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-g', '--galaxy', dest='gal', help='galaxy number')

options, args = parser.parse_args()

gal = str(args[0])


pp(gal)

