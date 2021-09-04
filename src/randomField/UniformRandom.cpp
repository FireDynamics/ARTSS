
#include "./UniformRandom.h"


Field UniformRandom::random_field(size_t size) {
    Field ret(size);

    for (size_t i=0; i < size; ++i) {
        ret[i] = m_dist(m_mt);
    }

    return ret;
}
