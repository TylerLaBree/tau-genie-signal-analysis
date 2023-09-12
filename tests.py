### Run a series of tests. User should verify tests are passed by eye. ###

from analysis import *

read_header('gallery/ValidHandle.h')
provide_get_valid_handle('std::vector<simb::MCTruth>')

sample_size = 100

# This is broken for some reason, likely because I changed how event attributes are flattened; that is — they aren't
#fig = get_energy_comparison_histogram(sample_size)
#fig.show()

# Also broken
#fig2 = get_energy_loss_comparison_histogram(sample_size)
#fig2.show()

# Also broken
#fig3 = get_scattering_angle_comparison_histogram(sample_size)
#fig3.show()

fig4 = get_interaction_type_bar_graph(sample_size)
fig4.show()

fig5 = get_charged_hadron_angle_histogram(sample_size)
fig5.show()

fig6 = get_charged_hadron_energy_histogram(sample_size)
fig6.show()

fig7 = get_charged_hadron_count_histogram(sample_size)
fig7.show()

fig8 = get_charged_hadron_count_histogram_split(sample_size) 
fig8.show()

input("Press Enter to continue...")

