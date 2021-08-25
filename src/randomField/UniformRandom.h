
#ifndef ARTSS_SOURCE_URANDOM_H_
#define ARTSS_SOURCE_URANDOM_H_


#include <random>
#include <iostream>

#include "../interfaces/IRandomField.h"


class UniformRandom: public IRandomField {
 public:
    UniformRandom(real range, real step_size) :
        m_seed(std::random_device()()),
        m_steps(range / step_size) {
        init_dist();
    }

    UniformRandom(real range, real step_size, int seed) :
        m_seed(seed), m_steps(range / step_size) {
        init_dist();
    }

    ~UniformRandom() = default;

    Field random_field(size_t size) override;

 private:
     void init_dist() {
         m_mt = std::mt19937(m_seed);
         m_dist = std::uniform_int_distribution<int>(-m_steps, m_steps);
     }

     const int m_seed;
     const real m_steps;
     std::mt19937 m_mt;
     std::uniform_int_distribution<int> m_dist;
};


#endif
