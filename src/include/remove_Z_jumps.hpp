#ifndef __REMOVE_Z_JUMPS__
#define __REMOVE_Z_JUMPS__


#include <neurostr/core/neuron.h>
#include <neurostr/core/node.h>

#include <neurostr/selector/selector.h>
#include <neurostr/selector/node_selector.h>
#include <neurostr/selector/neuron_selector.h>

void remove_Z_jumps_neurite(neurostr::Neurite& n);
void remove_Z_jumps_neuron(neurostr::Neuron& n);

#endif // __REMOVE_Z_JUMPS__

