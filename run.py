### Example run script ###

from datetime import datetime
from analysis import *

read_header('gallery/ValidHandle.h')
provide_get_valid_handle('std::vector<simb::MCTruth>')
today = datetime.today().strftime('%Y-%m-%d')

sample_size = 10000

fig = get_charged_hadron_count_histogram_split(sample_size) 
fig.show()
fig.savefig("/nashome/t/tlabree/img/"+today+"/charged_hadron_count.svg")

input("Press Enter to continue...")
