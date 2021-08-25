
#include "./UniformRandom.h"


Field UniformRandom::random_field() {
    Field ret(UNKNOWN_FIELD);
    size_t size = ret.get_size();

    for (size_t i=0; i < size; ++i) {
        ret[i] = m_dist(m_mt);
    }

    return ret;
}

